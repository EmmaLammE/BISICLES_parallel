import re
import sys

def extract_runtimes(file_path):
    runtimes = []
    with open(file_path, 'r') as file:
        for line in file:
            match = re.search(r'Total simulation time using PFASST: ! (\d+\.\d+) ! sec', line)
            if match:
                runtimes.append(float(match.group(1)))
    return runtimes

def calculate_statistics(runtimes):
    if not runtimes:
        return None, None, None
    mean_runtime = sum(runtimes) / len(runtimes)
    max_runtime = max(runtimes)
    min_runtime = min(runtimes)
    return mean_runtime, max_runtime, min_runtime

def main(file_path):
    runtimes = extract_runtimes('./sbatch_outputs/'+file_path)
    mean_runtime, max_runtime, min_runtime = calculate_statistics(runtimes)
    print("Mean Runtime: {} sec".format(mean_runtime))
    print("Max Runtime: {} sec".format(max_runtime))
    print("Min Runtime: {} sec".format(min_runtime))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_your_file.out>")
        sys.exit(1)
    file_path = sys.argv[1]
    main(file_path)
