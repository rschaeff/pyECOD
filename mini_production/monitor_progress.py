#!/usr/bin/env python3
"""
Monitor Mini PyECOD Production Progress

Real-time monitoring of job progress and results
Shows completion rates, identifies ready results for import

Usage:
    python monitor_progress.py                    # One-time status
    python monitor_progress.py --watch            # Continuous monitoring  
    python monitor_progress.py --ready-import     # Show ready for import
    python monitor_progress.py --stats            # Detailed statistics
"""

import os
import sys
import time
import subprocess
import argparse
from pathlib import Path
from typing import Dict, List, Tuple
from datetime import datetime
import yaml
import json

class MiniProgressMonitor:
    """Monitor mini PyECOD production progress"""
    
    def __init__(self, config_path: str = "config.local.yml"):
        self.config = self._load_config(config_path)
        
    def _load_config(self, config_path: str) -> Dict:
        """Load configuration"""
        try:
            with open(config_path, 'r') as f:
                return yaml.safe_load(f)
        except FileNotFoundError:
            # Use defaults if config not found
            return {
                "paths": {
                    "batch_base_dir": "/data/ecod/pdb_updates/batches",
                    "tracking_db": "/tmp/mini_production_status.db"
                }
            }
    
    def get_slurm_job_status(self) -> Dict[str, int]:
        """Get current SLURM job status"""
        
        try:
            # Get running jobs for current user
            result = subprocess.run(
                ["squeue", "-u", os.getenv("USER", "unknown"), "--noheader", "-o", "%i,%T"],
                capture_output=True, text=True
            )
            
            status_counts = {"RUNNING": 0, "PENDING": 0, "COMPLETING": 0}
            
            if result.stdout.strip():
                for line in result.stdout.strip().split('\n'):
                    if ',' in line:
                        job_id, status = line.strip().split(',')
                        if status in status_counts:
                            status_counts[status] += 1
                        elif status in ["CONFIGURING", "RESIZING"]:
                            status_counts["PENDING"] += 1
            
            return status_counts
            
        except Exception as e:
            print(f"Warning: Could not get SLURM status: {e}")
            return {"RUNNING": 0, "PENDING": 0, "COMPLETING": 0}
    
    def scan_batch_progress(self) -> Dict:
        """Scan all batches for progress information"""
        
        batch_base = Path(self.config["paths"]["batch_base_dir"])
        
        if not batch_base.exists():
            return {"error": f"Batch directory not found: {batch_base}"}
        
        batch_progress = {}
        total_proteins = 0
        total_completed = 0
        
        for batch_dir in batch_base.iterdir():
            if not batch_dir.is_dir():
                continue
            
            batch_name = batch_dir.name
            domains_dir = batch_dir / "domains"
            mini_domains_dir = batch_dir / "mini_domains"
            
            if not domains_dir.exists():
                continue
            
            # Count input files (proteins to process)
            input_files = list(domains_dir.glob("*.develop291.domain_summary.xml"))
            
            # Count output files (completed proteins)
            completed_files = []
            if mini_domains_dir.exists():
                completed_files = list(mini_domains_dir.glob("*.mini.domains.xml"))
            
            batch_stats = {
                "input_proteins": len(input_files),
                "completed_proteins": len(completed_files),
                "completion_rate": len(completed_files) / max(1, len(input_files)) * 100,
                "pending_proteins": len(input_files) - len(completed_files)
            }
            
            batch_progress[batch_name] = batch_stats
            total_proteins += len(input_files)
            total_completed += len(completed_files)
        
        # Overall statistics
        overall_stats = {
            "total_proteins": total_proteins,
            "total_completed": total_completed,
            "overall_completion_rate": total_completed / max(1, total_proteins) * 100,
            "total_pending": total_proteins - total_completed,
            "batches_scanned": len(batch_progress)
        }
        
        return {
            "overall": overall_stats,
            "batches": batch_progress,
            "timestamp": datetime.now().isoformat()
        }
    
    def get_ready_for_import(self) -> List[Dict]:
        """Get batches with completed results ready for import"""
        
        progress_data = self.scan_batch_progress()
        ready_batches = []
        
        for batch_name, stats in progress_data.get("batches", {}).items():
            if stats["completed_proteins"] > 0:
                ready_batches.append({
                    "batch_name": batch_name,
                    "completed_proteins": stats["completed_proteins"],
                    "completion_rate": stats["completion_rate"]
                })
        
        # Sort by number of completed proteins
        ready_batches.sort(key=lambda x: x["completed_proteins"], reverse=True)
        
        return ready_batches
    
    def print_status_summary(self):
        """Print a comprehensive status summary"""
        
        print("🔍 Mini PyECOD Production Status")
        print("=" * 60)
        print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        
        # SLURM job status
        slurm_status = self.get_slurm_job_status()
        print("📊 SLURM Jobs:")
        print(f"  Running:    {slurm_status['RUNNING']:>6}")
        print(f"  Pending:    {slurm_status['PENDING']:>6}")
        print(f"  Completing: {slurm_status['COMPLETING']:>6}")
        print(f"  Total:      {sum(slurm_status.values()):>6}")
        print()
        
        # Batch progress
        progress_data = self.scan_batch_progress()
        
        if "error" in progress_data:
            print(f"❌ Error: {progress_data['error']}")
            return
        
        overall = progress_data["overall"]
        print("🎯 Overall Progress:")
        print(f"  Total proteins:     {overall['total_proteins']:>8,}")
        print(f"  Completed:          {overall['total_completed']:>8,}")
        print(f"  Pending:            {overall['total_pending']:>8,}")
        print(f"  Completion rate:    {overall['overall_completion_rate']:>7.1f}%")
        print(f"  Batches scanned:    {overall['batches_scanned']:>8}")
        print()
        
        # Top completing batches
        batches = progress_data["batches"]
        top_batches = sorted(batches.items(), 
                           key=lambda x: x[1]["completed_proteins"], 
                           reverse=True)[:10]
        
        print("🏆 Top Completing Batches:")
        print(f"{'Batch Name':<40} {'Completed':<10} {'Total':<8} {'Rate':<8}")
        print("-" * 70)
        for batch_name, stats in top_batches:
            if stats["completed_proteins"] > 0:
                print(f"{batch_name:<40} {stats['completed_proteins']:<10} "
                      f"{stats['input_proteins']:<8} {stats['completion_rate']:<7.1f}%")
        print()
        
        # Import readiness
        ready_batches = self.get_ready_for_import()
        ready_total = sum(b["completed_proteins"] for b in ready_batches)
        
        print(f"✅ Ready for Import: {ready_total:,} results across {len(ready_batches)} batches")
        
    def print_detailed_stats(self):
        """Print detailed statistics"""
        
        progress_data = self.scan_batch_progress()
        
        if "error" in progress_data:
            print(f"❌ Error: {progress_data['error']}")
            return
        
        print("📈 Detailed Batch Statistics")
        print("=" * 80)
        
        batches = progress_data["batches"]
        
        # Group by completion status
        not_started = []
        in_progress = []
        high_completion = []
        
        for batch_name, stats in batches.items():
            completion_rate = stats["completion_rate"]
            if completion_rate == 0:
                not_started.append((batch_name, stats))
            elif completion_rate < 80:
                in_progress.append((batch_name, stats))
            else:
                high_completion.append((batch_name, stats))
        
        # Print categories
        print(f"\n🔴 Not Started ({len(not_started)} batches):")
        for batch_name, stats in sorted(not_started, key=lambda x: x[1]["input_proteins"], reverse=True)[:5]:
            print(f"  {batch_name:<40} {stats['input_proteins']:>8} proteins")
        
        print(f"\n🟡 In Progress ({len(in_progress)} batches):")
        for batch_name, stats in sorted(in_progress, key=lambda x: x[1]["completion_rate"], reverse=True):
            print(f"  {batch_name:<40} {stats['completed_proteins']:>6}/{stats['input_proteins']:<6} "
                  f"({stats['completion_rate']:>5.1f}%)")
        
        print(f"\n🟢 High Completion ({len(high_completion)} batches):")
        for batch_name, stats in sorted(high_completion, key=lambda x: x[1]["completion_rate"], reverse=True):
            print(f"  {batch_name:<40} {stats['completed_proteins']:>6}/{stats['input_proteins']:<6} "
                  f"({stats['completion_rate']:>5.1f}%)")
    
    def print_import_readiness(self):
        """Print information about results ready for import"""
        
        ready_batches = self.get_ready_for_import()
        
        print("📥 Import Readiness Report")
        print("=" * 60)
        
        if not ready_batches:
            print("No completed results found yet.")
            return
        
        total_ready = sum(b["completed_proteins"] for b in ready_batches)
        print(f"Total ready for import: {total_ready:,} proteins\n")
        
        print(f"{'Batch Name':<40} {'Ready':<8} {'Rate':<8}")
        print("-" * 60)
        
        for batch in ready_batches:
            print(f"{batch['batch_name']:<40} {batch['completed_proteins']:<8} "
                  f"{batch['completion_rate']:<7.1f}%")
        
        print(f"\nSuggested import commands:")
        print(f"# Import all ready results:")
        print(f"python import_results.py --import-all")
        print(f"")
        print(f"# Import specific high-completion batch:")
        if ready_batches:
            top_batch = ready_batches[0]
            print(f"python import_results.py --batch-name {top_batch['batch_name']}")
        print(f"")
        print(f"# Import limited number for testing:")
        print(f"python import_results.py --import-all --limit 100")
    
    def watch_progress(self, interval: int = 60):
        """Watch progress continuously"""
        
        print("👀 Watching Mini PyECOD Progress (Ctrl+C to stop)")
        print(f"Update interval: {interval} seconds\n")
        
        try:
            while True:
                # Clear screen
                os.system('clear' if os.name == 'posix' else 'cls')
                
                self.print_status_summary()
                
                print(f"\nNext update in {interval} seconds...")
                time.sleep(interval)
                
        except KeyboardInterrupt:
            print("\n👋 Monitoring stopped")


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description='Monitor Mini PyECOD Production Progress'
    )
    
    parser.add_argument('--watch', action='store_true',
                       help='Watch progress continuously')
    parser.add_argument('--interval', type=int, default=60,
                       help='Update interval for watch mode (seconds)')
    parser.add_argument('--ready-import', action='store_true',
                       help='Show batches ready for import')
    parser.add_argument('--stats', action='store_true',
                       help='Show detailed statistics')
    parser.add_argument('--config', type=str, default='config.local.yml',
                       help='Config file path')
    
    args = parser.parse_args()
    
    # Initialize monitor
    monitor = MiniProgressMonitor(args.config)
    
    if args.watch:
        monitor.watch_progress(args.interval)
    elif args.ready_import:
        monitor.print_import_readiness()
    elif args.stats:
        monitor.print_detailed_stats()
    else:
        monitor.print_status_summary()


if __name__ == "__main__":
    main()
