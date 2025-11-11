#!/usr/bin/env python
import glafic
import numpy as np
import psutil
import shutil
import os
import time
import requests
import json
import pandas as pd
import re
import sys
import csv
import multiprocessing
from itertools import product
import uuid


# ==== Config ====
m = [round(x, 5) for x in np.linspace(0.05, 0.1, 2)]
m_lens = 2
m_param = 5
n = [round(x, 5) for x in np.linspace(0, 360, 100)]
n_lens = 2
n_param = 6

# --- File Paths ---
model_output_dir = '/Volumes/T7 Shield/Simulations/Output' # Model Save Path
base_results_path = '/Volumes/T7 Shield/Simulations/Input' # Path with all input files
log_file_path = '/Volumes/T7 Shield/Simulations/Run/simulation_log.txt' # Log file path
restart_file_path = os.path.join(os.path.dirname(log_file_path), 'simulation_restart_state_var.json') # Restart file
obs_point_file = base_results_path + '/pos+flux_point.dat'  # Observation file path ( Should include all available data )
constraint_file = base_results_path + '/pos_point.dat'  # Constraint file path ( Same as obs_point but changed for the specific constraint )
prior_file = None # Path to prior file
input_py_file = '/Volumes/T7 Shield/Simulations/Input/input.py' # Input Python file path
sim_name = 'Sim 1' # Name of Sim
model = 'SIE' # Model Type

# --- Performance & Resource Management ---
NUM_PROCESSORS = 6 # Number of CPU cores to use
CHUNK_SIZE = 100000  # Sims split into chuncks of this size
CONSOLE_UPDATE_INTERVAL = 100 # Update console output every N iterations

# Ensure the output directories exist
os.makedirs(model_output_dir, exist_ok=True)
os.makedirs(base_results_path, exist_ok=True)


# ==== Helper Functions ====
def get_cpu_usage():
    return psutil.cpu_percent(interval=None)

def get_memory_usage():
    return psutil.virtual_memory().percent

def get_dir_size(directory):
    total_size = 0
    try:
        for dirpath, dirnames, filenames in os.walk(directory):
            for f in filenames:
                fp = os.path.join(dirpath, f)
                if not os.path.islink(fp): total_size += os.path.getsize(fp)
    except FileNotFoundError: return 0
    return total_size

def get_disk_usage_info(directory, avg_model_size, models_remaining):
    try:
        disk_root = os.path.abspath(directory)
        while not os.path.ismount(disk_root):
            parent = os.path.dirname(disk_root)
            if parent == disk_root: break
            disk_root = parent
        total, used, free = shutil.disk_usage(disk_root)
        used_percent = (used / total) * 100
        free_gb = free / (1024**3)
        projected_usage = avg_model_size * models_remaining
        if projected_usage > free: threat = "HIGH"
        elif projected_usage > 0.8 * free: threat = "MEDIUM"
        else: threat = "LOW"
        return {"used_percent": used_percent, "free_gb": free_gb, "threat_level": threat}
    except (FileNotFoundError, Exception):
        return {"used_percent": 0, "free_gb": 0, "threat_level": "Unavailable"}

def save_restart_state(path, i, j, chunk_number):
    state = {'i': i, 'j': j, 'chunk_number': chunk_number}
    with open(path, 'w') as f:
        json.dump(state, f, indent=4)
    with open(sys.__stdout__.fileno(), 'w', closefd=False) as real_stdout:
        real_stdout.write(f"\n✅ Simulation state saved to {path} at indices (i={i}, j={j}), chunk={chunk_number}.\n")

def load_restart_state(path):
    if not os.path.exists(path):
        with open(sys.__stdout__.fileno(), 'w', closefd=False) as real_stdout:
            real_stdout.write("ℹ️ Restart file not found. Starting a fresh simulation.\n")
        return 0, 0, 1
    try:
        with open(path, 'r') as f:
            state = json.load(f)
            i, j = state.get('i', 0), state.get('j', 0)
            chunk = state.get('chunk_number', 1)
            with open(sys.__stdout__.fileno(), 'w', closefd=False) as real_stdout:
                real_stdout.write(f"✅ Restart file found. Resuming from (i={i}, j={j}), chunk={chunk}.\n")
            return i, j, chunk
    except (json.JSONDecodeError, KeyError):
        with open(sys.__stdout__.fileno(), 'w', closefd=False) as real_stdout:
            real_stdout.write(f"⚠️ Could not read restart file. Starting fresh.\n")
        return 0, 0, 1

def get_csv_filename(chunk_number):
    # sim_name = model_output_dir.strip('/').split('/')[-1]
    return os.path.join(model_output_dir, f"{sim_name}_summary.csv" if chunk_number == 1 else f"{sim_name}_summary_{chunk_number}.csv")

def write_to_csv(df, chunk_number):
    csv_file = get_csv_filename(chunk_number)
    header = not os.path.exists(csv_file)
    df.to_csv(csv_file, mode='a', index=False, header=header)

def calculate_chunk_number(iteration_count):
    return ((iteration_count - 1) // CHUNK_SIZE) + 1

def rms_extract(model_ver, model_path, temp_input_py_file):
    global pos_rms, mag_rms, chi2_value
    
    # --- CORRECTED: Use os.path.join for robust path construction ---
    opt_result_file = os.path.join(model_path, f'{model_ver}_optresult.dat')
    
    # Load the data
    with open(opt_result_file, 'r') as file:
        opt_result = file.readlines()

    # Find the last line with 'optimize' in it
    last_optimize_index = None
    for idx in range(len(opt_result) - 1, -1, -1):
        if 'optimize' in opt_result[idx]:
            last_optimize_index = idx
            last_optimize_line = opt_result[idx]
            break
    if last_optimize_index is None:
        raise ValueError("No line with 'optimize' found in the file.")

    # Extract everything after the last 'optimize' line
    opt_result = opt_result[last_optimize_index + 1:]

    # Count the number of lines that start with 'lens'
    lens_count = sum(1 for line in opt_result if line.startswith('lens'))
    lens_params_dict = {}
    lens_params = []
    for line in opt_result:
        if line.startswith('lens'):
            parts = re.split(r'\s+', line.strip())
            lens_name = parts[1]
            params = [float(x) for x in parts[2:]]
            lens_params_dict[lens_name] = params
            lens_params.append((lens_name, params))
    if lens_params:
        for i in range(len(lens_params)):
            lens_name, params = lens_params[i]
            lens_params_dict[lens_name] = params[1:]

    source_params = []
    for line in opt_result:
        if line.startswith('point'):
            parts = re.split(r'\s+', line.strip())
            params = [float(x) for x in parts[1:]]
            source_params.append(params)
    
    chi2_line = next((line for line in opt_result if 'chi^2' in line), None)
    if chi2_line is None:
        raise ValueError("No line with 'chi2' found in the file.")
    chi2_value = float(chi2_line.split('=')[-1].strip().split()[0])
    num_lens_profiles = len(lens_params_dict)

    df = pd.DataFrame()
    rows = []
    max_param_len = 0
    for lens_name, params in lens_params_dict.items():
        row = {'Lens Name': lens_name}
        for i, val in enumerate(params): row[f'param{i+1}'] = val
        rows.append(row)
        if len(params) > max_param_len: max_param_len = len(params)
    columns = ['Lens Name'] + [f'param{i+1}' for i in range(max_param_len)]
    df = pd.DataFrame(rows, columns=columns)

    with open(temp_input_py_file, 'r') as file:
        py = file.readlines()

    set_lens_lines = [line for line in py if line.startswith('glafic.set_lens(')]
    if not set_lens_lines: raise ValueError("No lines starting with 'glafic.set_lens(' found in the file.")
    set_lens_params = []
    for line in set_lens_lines:
        match = re.search(r'set_lens\((.*?)\)', line)
        if match:
            params_str = match.group(1)
            params = [param.strip() for param in params_str.split(',')]
            set_lens_params.append(params)
        else: raise ValueError(f"No valid parameters found in line: {line.strip()}")

    set_lens_dict = {}
    for params in set_lens_params:
        if len(params) < 3: raise ValueError(f"Not enough parameters found in line: {params}")
        lens_name = params[1].strip("'\"")
        lens_params = [float(x) for x in params[2:]]
        set_lens_dict[lens_name] = lens_params
    if set_lens_dict:
        for lens_name, params in set_lens_dict.items():
            set_lens_dict[lens_name] = params[1:]

    df_input = pd.DataFrame()
    rows_input = []
    max_param_len_input = 0
    for lens_name, params in set_lens_dict.items():
        row = {'Lens Name': lens_name}
        for i, val in enumerate(params): row[f'param{i+1}'] = val
        rows_input.append(row)
        if len(params) > max_param_len_input: max_param_len_input = len(params)
    columns_input = ['Lens Name'] + [f'param{i+1}' for i in range(max_param_len_input)]
    df_input = pd.DataFrame(rows_input, columns=columns_input)
    
    set_flag_lines = [line for line in py if line.startswith('glafic.setopt_lens(')]
    if not set_flag_lines: raise ValueError("No lines starting with 'glafic.setopt_lens(' found in the file.")
    set_flag_params = []
    for line in set_flag_lines:
        match = re.search(r'setopt_lens\((.*?)\)', line)
        if match:
            params_str = match.group(1)
            params = [param.strip() for param in params_str.split(',')]
            set_flag_params.append(params)
        else: raise ValueError(f"No valid parameters found in line: {line.strip()}")
    
    set_flag_dict = {}
    for params in set_flag_params:
        if len(params) < 2: raise ValueError(f"Not enough parameters found in line: {params}")
        lens_index = params[0].strip("'\"")
        lens_name = None
        for lens_params in set_lens_params:
            if lens_params[0].strip("'\"") == lens_index:
                lens_name = lens_params[1].strip("'\"")
                break
        if lens_name is None: raise ValueError(f"Lens name for index {lens_index} not found")
        flag = ','.join(params[1:])
        set_flag_dict[lens_name] = flag
    if set_flag_dict:
        for lens_name, flag in set_flag_dict.items():
            flag_parts = flag.split(',')
            set_flag_dict[lens_name] = ','.join(flag_parts[1:])
    
    df_flag = pd.DataFrame()
    rows_flag = []
    max_flag_len = 0
    for flag in set_flag_dict.values():
        flag_parts = flag.split(',')
        if len(flag_parts) > max_flag_len: max_flag_len = len(flag_parts)
    for lens_name, flag in set_flag_dict.items():
        flag_parts = flag.split(',')
        row = {'Lens Name': lens_name}
        for i, val in enumerate(flag_parts): row[f'flag{i+1}'] = val
        rows_flag.append(row)
    columns_flag = ['Lens Name'] + [f'flag{i+1}' for i in range(max_flag_len)]  
    df_flag = pd.DataFrame(rows_flag, columns=columns_flag)
    
    dfs = []
    for i in range(num_lens_profiles):
        lens_name = df['Lens Name'][i]
        model_type = None
        for m in model_list:
            if m.lower() == lens_name.lower():
                model_type = m
                break
        if model_type is None: continue
        symbols = model_params[model_type][:7]
        row_input = pd.DataFrame([df_input.iloc[i, 1:8].values], columns=symbols)
        row_output = pd.DataFrame([df.iloc[i, 1:8].values], columns=symbols)
        row_flags = pd.DataFrame([df_flag.iloc[i, 1:8].values], columns=symbols)
        lens_df = pd.concat([row_input.assign(Type='Input'), row_output.assign(Type='Output'), row_flags.assign(Type='Flag')], ignore_index=True)
        lens_df.insert(0, 'Lens Name', lens_name)
        cols = lens_df.columns.tolist()
        cols.insert(1, cols.pop(cols.index('Type')))
        lens_df = lens_df[cols]
        dfs.append(lens_df)
    
    columnn_names = ['x', 'y', 'mag', 'pos_err', 'mag_err', '1', '2', '3']
    obs_point = pd.read_csv(obs_point_file, delim_whitespace=True, header=None, skiprows=1, names=columnn_names)
    
    # --- CORRECTED: Also use os.path.join for the point file path ---
    out_point_file = os.path.join(model_path, f'{model_ver}_point.dat')
    out_point = pd.read_csv(out_point_file, delim_whitespace=True, header=None, skiprows=1, names=columnn_names)
    
    out_point.drop(columns=['mag_err', '1', '2', '3'], inplace=True)
    mask = abs(out_point['mag']) >= 1
    out_point = out_point[mask[:len(out_point)]].reset_index(drop=True)
    out_point['x_diff'] = abs(out_point['x'] - obs_point['x'])
    out_point['y_diff'] = abs(out_point['y'] - obs_point['y'])
    out_point['mag_diff'] = abs(abs(out_point['mag']) - abs(obs_point['mag']))
    out_point['pos_sq'] = np.sqrt((out_point['x_diff']**2 + out_point['y_diff']**2).astype(float))
    pos_rms = np.average(out_point['pos_sq'])
    mag_rms = np.average(np.sqrt((out_point['mag_diff']**2).astype(float)))
    return pos_rms, mag_rms, dfs, chi2_value, source_params

pow_params = ['$z_{s,fid}$', 'x', 'y', 'e', '$θ_{e}$', '$r_{Ein}$', '$\gamma$ (PWI)']
sie_params = ['$\sigma$', 'x', 'y', 'e', '$θ_{e}$', '$r_{core}$', 'NaN']
nfw_params = ['M', 'x', 'y', 'e', '$θ_{e}$', 'c or $r_{s}$', 'NaN']
ein_params = ['M', 'x', 'y', 'e', '$θ_{e}$', 'c or $r_{s}$', r'$\alpha_{e}$']
shear_params = ['$z_{s,fid}$', 'x', 'y', '$\gamma$', '$θ_{\gamma}$', 'NaN', '$\kappa$']
sersic_params = ['$M_{tot}$', 'x', 'y', 'e', '$θ_{e}$', '$r_{e}$', '$n$']
cored_sie_params = ['M', 'x', 'y', 'e', '$θ_{e}$', '$r_{core}$', 'NaN']
mpole_params = ['$z_{s,fid}$', 'x', 'y', '$\epsilon$', '$θ_{m}$', 'm', 'n']
model_list = ['POW', 'SIE', 'ANFW', 'EIN', 'PERT', 'SERS', 'MPOLE']
model_params = {'POW': pow_params, 'SIE': sie_params, 'ANFW': nfw_params, 'EIN': ein_params, 'PERT': shear_params, 'SERS': sersic_params, 'MPOLE': mpole_params}


# ==== Worker Function for Multiprocessing ====
def run_single_model(params):

    i, j, m_val, n_val = params  # Unpack 2 variables and their indices

    model_name = model + f'_{m_val}_{n_val}'
    
    unique_id = uuid.uuid4()
    temp_input_py_file = os.path.join(os.path.dirname(input_py_file), f"temp_input_{unique_id}.py")
    
    try:
        shutil.copy(input_py_file, temp_input_py_file)

        with open(temp_input_py_file, 'r') as f:
            lines = f.readlines()

        for line_num, line in enumerate(lines):
            if line.strip().startswith('path ='):
                lines[line_num] = f"path = '{model_output_dir}/{model_name}'\n"
            elif line.strip().startswith('constraint_file ='):
                lines[line_num] = f"constraint_file = '{constraint_file}'\n"
        
        lens_lines_indices = [idx for idx, line in enumerate(lines) if line.strip().startswith('glafic.set_lens(')]
        if not lens_lines_indices:
            raise ValueError(f"CRITICAL: No 'glafic.set_lens' lines found in '{temp_input_py_file}'.")

        params_to_set = [(m_lens, m_param, m_val), (n_lens, n_param, n_val)]

        for lens_target, param_index, value in params_to_set:
            if lens_target > len(lens_lines_indices):
                raise IndexError(f"Config Error: m/n_lens is {lens_target}, but only {len(lens_lines_indices)} 'glafic.set_lens' lines found.")
            
            line_index_in_file = lens_lines_indices[lens_target - 1]
            line = lines[line_index_in_file]
            match = re.search(r'\((.*)\)', line)
            if not match: raise ValueError(f"Could not parse parameters from line: {line}")
            
            parts = [p.strip() for p in match.group(1).split(',')]
            target_index_in_parts = 2 + (param_index - 1)

            if target_index_in_parts >= len(parts):
                raise IndexError(f"Config Error: m/n_param is {param_index}, which is out of range for the parameters in line:\n>>> {line.strip()}")
            
            parts[target_index_in_parts] = str(value)
            lines[line_index_in_file] = re.sub(r'\(.*\)', f'({", ".join(parts)})', line)

        with open(temp_input_py_file, 'w') as f:
            f.writelines(lines)
    
        os.system(f'python3 "{temp_input_py_file}"')

        pos_rms, mag_rms, dfs, chi2, source = rms_extract(model_name, model_output_dir, temp_input_py_file)
        if pos_rms == -1: raise IOError("rms_extract failed.")
        
        num_images = 0
        out_point_file = os.path.join(model_output_dir, f'{model_name}_point.dat')
        if os.path.exists(out_point_file):
            with open(out_point_file, 'r') as f: num_images = sum(1 for line in f if line.strip()) - 1
            
        result_dict = {'m': m_val, 'n': n_val, 'num_images': num_images, 'pos_rms': pos_rms, 'mag_rms': mag_rms, 'chi2': chi2}
        

        num_lenses = model_name.count('_') - 1
        macro_model_params = model_name.strip().split('_')[0]
        macro_columns = model_params[macro_model_params]

        if num_lenses > 1:
            for i in range(1, num_lenses):
                micro_model_params = model_name.strip().split('_')[i]
                if micro_model_params == 'SHEAR': micro_model_params = 'PERT'
                micro_columns = model_params[micro_model_params]
                micro_columns = [f'{col}_{i}' for col in micro_columns]

        macro_cols = list(dict.fromkeys(macro_columns[:7])) if isinstance(macro_columns, list) else [macro_columns]
        if num_lenses > 1:
            micro_cols = list(dict.fromkeys(micro_columns[:7])) if isinstance(micro_columns, list) else [micro_columns]

        def _safe_get(dfs_idx, col):
            if len(dfs) > dfs_idx and col in dfs[dfs_idx].columns:
                return dfs[dfs_idx][col].iloc[1]
            return 0

        for col in macro_cols: result_dict[col] = _safe_get(0, col)
        if num_lenses > 1:
            for len_num in range(1, num_lenses):
                for col in micro_cols:
                    stripped_col = col.strip(f'_{len_num}')
                    result_dict[col] = _safe_get(1, stripped_col)
        
        result_dict['source_x'] = source[0][1] if source and len(source) > 0 else 0
        result_dict['source_y'] = source[0][2] if source and len(source) > 0 else 0

    except Exception as e:
        import traceback
        error_info = f"Error in model {model_name}: {e}\n{traceback.format_exc()}"
        result_dict = {'error': error_info}
    finally:
        if os.path.exists(temp_input_py_file):
            os.remove(temp_input_py_file)
            
        for suffix in ['_point.dat', '_optresult.dat']:
            file_to_delete = os.path.join(model_output_dir, f"{model_name}{suffix}")
            if os.path.exists(file_to_delete):
                os.remove(file_to_delete)
            
    return (i, j, result_dict)


# ==== Main Process ====
def main():
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__

    print(f"Starting 2-variable simulation with {NUM_PROCESSORS} processes.")
    all_tasks = [(i, j, m[i], n[j]) for i, j in product(range(len(m)), range(len(n)))]
    total_iterations = len(all_tasks)
    print(f"Total iterations to run: {total_iterations}")

    start_i, start_j, chunk_number = load_restart_state(restart_file_path)
    start_index = 0
    if start_i > 0 or start_j > 0:
        try:
            last_idx = [task[:2] for task in all_tasks].index((start_i, start_j))
            start_index = last_idx + 1
        except ValueError:
            print("Warning: Restart indices not found. Starting from scratch.")
    tasks_to_run = all_tasks[start_index:]
    iterations_done = start_index
    
    if not tasks_to_run:
        print("🎉 Simulation already completed!")
        if os.path.exists(restart_file_path): os.remove(restart_file_path)
        return
        
    print(f"Resuming from iteration {iterations_done}. {len(tasks_to_run)} tasks remaining.")

    start_time = time.time()
    last_state_saved_at_iteration = iterations_done
    last_i, last_j = start_i, start_j
    avg_model_size, initial_dir_size = 0, get_dir_size(base_results_path)
    disk_info = {"used_percent": 0, "free_gb": 0, "threat_level": "Calculating..."}

    try:
        with multiprocessing.Pool(processes=NUM_PROCESSORS) as pool:
            results_iterator = pool.imap_unordered(run_single_model, tasks_to_run)
            
            for i_res, j_res, result_dict in results_iterator:
                if 'error' in result_dict:
                    print(f"Warning: A worker failed with error: {result_dict['error']}")
                    continue
                
                last_i, last_j = i_res, j_res
                current_chunk = calculate_chunk_number(iterations_done + 1)
                write_to_csv(pd.DataFrame([result_dict]), current_chunk)
                iterations_done += 1

                if iterations_done % CONSOLE_UPDATE_INTERVAL == 0 or iterations_done == total_iterations:
                    percentage_complete = (iterations_done / total_iterations) * 100
                    print(f"Progress: {iterations_done} / {total_iterations} ({percentage_complete:.2f}%)")

                models_this_run = iterations_done - start_index
                
                if models_this_run > 0 and models_this_run % 100 == 0:
                    current_dir_size = get_dir_size(base_results_path)
                    avg_model_size = (current_dir_size - initial_dir_size) / models_this_run
                    disk_info = get_disk_usage_info(base_results_path, avg_model_size, total_iterations - iterations_done)
                
                percentage_complete = (iterations_done / total_iterations) * 100
                elapsed_time = time.time() - start_time
                avg_time_per_iteration = elapsed_time / models_this_run if models_this_run > 0 else 0
                approx_time_remaining = (total_iterations - iterations_done) * avg_time_per_iteration

                progress_info = {
                    'current_iteration': iterations_done,
                    'total_iterations': total_iterations,
                    'percentage_complete': percentage_complete,
                    'avg_time_per_iteration': avg_time_per_iteration,
                    'approx_time_remaining': approx_time_remaining,
                    'ram_usage_percent': get_memory_usage(),
                    'disk_threat_level': disk_info['threat_level'],
                    'cpu_usage_percent': get_cpu_usage(),
                    'disk_used_percent': disk_info['used_percent'],
                    'disk_free': round(disk_info['free_gb'], 2) if isinstance(disk_info['free_gb'], (int, float)) else disk_info['free_gb'],
                    'current_chunk': current_chunk
                }

                with open(log_file_path, 'w') as log_file:
                    log_file.write(f"Iteration {iterations_done}/{total_iterations}: {progress_info}\n")

                if iterations_done - last_state_saved_at_iteration >= 500:
                    save_restart_state(restart_file_path, last_i, last_j, current_chunk)
                    last_state_saved_at_iteration = iterations_done
    
    except KeyboardInterrupt:
        print("\n⚠️ Simulation interrupted by user.")
    except Exception as e:
        print(f"\n❌ An unexpected error occurred: {e}")
        import traceback
        traceback.print_exc()
    finally:
        print("\nFinalizing simulation state...")
        if iterations_done == total_iterations:
            print("🎉 Simulation completed successfully!")
            if os.path.exists(restart_file_path):
                os.remove(restart_file_path)
                print("🗑️ Removed restart file.")
        else:
            save_restart_state(restart_file_path, last_i, last_j, calculate_chunk_number(iterations_done))
        print("✅ Final state saved.")

if __name__ == "__main__":
    multiprocessing.set_start_method("spawn", force=True)
    main()