# Purpose: Run all MSA software on given datasets and capture the runtime and accuracy metrics for all software and datasets
# Author: Daniel Yuan
# Usage:
#     python run_benchmark.py
#     - Specify Parameters inside run_benchmark.py for ease of use
#     - Directories specified in parameters must exist before script is run

#Libraries
import os
import sys
import json
from subprocess import call
from timeit import timeit
import collections
import pdb

#Parameters:
source_dir = "./source/"
unaligned_dir = "./unaligned/"
software_dir = "./software/"
results_dir = "./results/"
#Helper Functions

# Unaligns all aligned files in the given source directory and places it in the target directory
# Returns array of paths to unaligned sequences
def unalign(source, target):
  unaligned_files = []
  source_files = sorted(os.listdir(source))
  for source_file in source_files:
    target_path = os.path.join(target, source_file) 
    source_path = os.path.join(source, source_file)
    if ( not os.path.isfile(target_path)):
      call("sed 's/-//g' " + source_path + " > " + target_path, shell=True)
    unaligned_files.append({"name": source_file, "source_path": source_path, "unalign_path": target_path})
  return unaligned_files

# Runs the given command on the shell and returns amount of time taken for the given command
def time_func(cmd):
  time_taken = timeit(stmt = "call('" + cmd + "', shell=True)", setup= "from subprocess import call", number = 1)
  return time_taken

# Makes a directory if it doesn't already exist
def make_dir(path):
  if (not os.path.exists(path)):
    os.makedirs(path)
  return path

# Runs FastSP and extracts results for a reference and 
def benchmark(ref, est, storage):
  benchmark_file = os.path.join(storage, "fastSP_results")
  jar_file = os.path.join(software_dir, "FastSP", "FastSP.jar")
  call("java -jar " + jar_file + " -r " + ref + " -e " + est + " > " + benchmark_file, shell=True)
  results = parse_benchmark(benchmark_file)
  return results

# Parse FastSP_results file to extract accuracy parameters
def parse_benchmark(result):
  with open(result) as f:
    lines = f.readlines()
  
  results = {}
  results["SP_score"] = lines[0].split(" ")[-1].split("\n")[0]
  results["Modeler"] = lines[1].split(" ")[-1].split("\n")[0]
  results["SPFN"] = lines[2].split(" ")[-1].split("\n")[0]
  results["SPFP"] = lines[3].split(" ")[-1].split("\n")[0]
  results["TC"] = lines[5].split(" ")[-1].split("\n")[0]
  return results

# Averages results for all the replicate datasets
def avg_results(storage):
  for software in storage.keys():
    cur_soft = storage[software]

    # Check if software storage is empty
    if not cur_soft:
      continue

    avg_time = 0;
    avg_SP = 0;
    avg_Modeler = 0;
    avg_SPFN = 0;
    avg_SPFP = 0;
    avg_TC = 0;
    for data in cur_soft.keys():
      results = cur_soft[data]['results']
      time = cur_soft[data]['time']
      avg_time += float(time)
      avg_SP += float(results['SP_score'])
      avg_Modeler += float(results['Modeler'])
      avg_SPFN += float(results['SPFN'])
      avg_SPFP += float(results['SPFP'])
      avg_TC += float(results['TC'])

    num_datasets = len(cur_soft.keys())
    avg_results = {
      "avg_time": float(avg_time)/num_datasets,
      "avg_SP": float(avg_SP)/num_datasets,
      "avg_Modeler": float(avg_Modeler)/num_datasets,
      "avg_SPFN": float(avg_SPFN)/num_datasets,
      "avg_SPFP": float(avg_SPFP)/num_datasets,
      "avg_TC": float(avg_TC)/num_datasets,
    }
    cur_soft["avg"] = avg_results 

# Saves any python object to json file
def save_json(storage, file_path):
  with open(file_path, 'w') as f:
    json.dump(storage, f)

#Software Functions: define functions for all software
def mafft_prog(ref, unalign):
  storage_dir = make_dir(os.path.join(results_dir, "mafft_proj"))
  exe = os.path.join(software_dir, "mafft", "scripts", "mafft")
  est = os.path.join(storage_dir, "estimated_alignment")
  time_taken = time_func(exe + " --retree 1 " + unalign + " > " + est)
  results = benchmark(ref, est, storage_dir)
  return {"time": time_taken, "results": results}

def mafft_iter(ref, unalign):
  storage_dir = make_dir(os.path.join(results_dir, "mafft_iter"))
  exe = os.path.join(software_dir, "mafft", "scripts", "mafft")
  est = os.path.join(storage_dir, "estimated_alignment")
  time_taken = time_func(exe + " --localpair --maxiterate 1000 " + unalign + " > " + est)
  results = benchmark(ref, est, storage_dir)
  return {"time": time_taken, "results": results}

def muscle(ref, unalign):
  storage_dir = make_dir(os.path.join(results_dir, "muscle"))
  exe = os.path.join(software_dir, "MUSCLE", "muscle")
  est = os.path.join(storage_dir, "estimated_alignment")
  time_taken = time_func(exe + " -in " + unalign + " -out " + est)
  results = benchmark(ref, est, storage_dir)
  return {"time": time_taken, "results": results}

def PRRN(ref, unalign):
  storage_dir = make_dir(os.path.join(results_dir, "PRRN"))
  exe = os.path.join(software_dir, "PRRN", "bin", "prrn")
  est = os.path.join(storage_dir, "estimated_alignment")
  time_taken = time_func(exe + " -FN " + unalign + " > " + est)
  results = benchmark(ref, est, storage_dir)
  return {"time": time_taken, "results": results}

def T_coffee(ref, unalign):
  storage_dir = make_dir(os.path.join(results_dir, "T-Coffee"))
  exe = os.path.join(software_dir, "T-coffee", "bin", "t_coffee")
  est = os.path.join(storage_dir, "estimated_alignment")
  time_taken = time_func(exe + " -in " + unalign + " -output fasta_aln -outfile " + est)
  results = benchmark(ref, est, storage_dir)
  return {"time": time_taken, "results": results}
  
#Main Function
def main():
  # Unalign all sequence in source folder and put them in unaligned folder
  unaligned_seq = unalign(source_dir, unaligned_dir)

  storage = {
    "mafft_prog": {},
    "mafft_iter": {},
    "PRRN": {},
    "T-COFFEE": {},
    "PSAlign": {},
    "MUSCLE": {}
  }
  # Run All software on all unaligned datasets
  for seq in unaligned_seq:
    ref = seq["source_path"]
    unalign_path = seq["unalign_path"]
    name = seq["name"]

    pdb.set_trace()

    #Add Results for all software for this sequence
    storage["mafft_prog"][name] = mafft_prog(ref, unalign_path)
 

  # Average Results for all datasets
  avg_results(storage)
  # Save results to json file for each viewing
  output_file = os.path.join(results_dir, "Benchmark_results.json")
  save_json(storage, output_file)

if __name__ == "__main__":
  main()
