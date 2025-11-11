#!/usr/bin/env python
import os
import subprocess
import numpy as np
import pandas as pd
import re
import sys
import shutil
import uuid
import time

# --- File Paths & Parameters ---
# These paths point to the SHARED filesystem for final results and initial inputs
model_output_base = '/home/rommulus/Projects/itng_lensing/Simulations/Output'
base_results_path = '/home/rommulus/Projects/itng_lensing/Simulations/Input/System_2'
obs_point_file = os.path.join(base_results_path, 'pos_point.dat')

# --- NEW: Get the fast LOCAL scratch directory path from the environment variable ---
# Fallback to /tmp if the variable isn't set for some reason (e.g., local testing)
SCRATCH_DIR = os.getenv('SCRATCH_DIR', '/tmp')

# --- Simulation Parameters ---
sim_name = 'Sim 8'
model = 'POW'
m = [round(x, 5) for x in np.linspace(1.5, 2.2, 100)]
n = [round(x, 5) for x in np.linspace(0.0, 0.9, 100)]
o = [round(x, 5) for x in np.linspace(0, 360, 100)]
m_lens, m_param = 1, 8
n_lens, n_param = 1, 5
o_lens, o_param = 1, 6
constraint_file = os.path.join(base_results_path, 'pos_point.dat')
prior_file = os.path.join(base_results_path, 'prior.dat')
input_py_file = os.path.join(base_results_path, 'pow_input.py')
time_delay, h0 = False, False

# --- Directory setup for final output on the shared filesystem ---
final_results_dir = os.path.join(model_output_base, 'individual_results')
progress_status_dir = os.path.join(model_output_base, 'progress_status')
os.makedirs(final_results_dir, exist_ok=True)
os.makedirs(progress_status_dir, exist_ok=True)

PROGRESS_UPDATE_INTERVAL = 100

# --- Parameter definitions ---
pow_params = ['$z_{s,fid}$', 'x', 'y', 'e', '$θ_{e}$', '$r_{Ein}$', '$\gamma$ (PWI)']
sie_params = ['$\sigma$', 'x', 'y', 'e', '$θ_{e}$', '$r_{core}$', 'NaN']
nfw_params = ['M', 'x', 'y', 'e', '$θ_{e}$', 'c or $r_{s}$', 'NaN']
pert_params = ['$z_{s,fid}$', 'x', 'y', '$\gamma$', '$θ_{\gamma}$', 'NaN', '$\kappa$']
model_params = {'POW': pow_params, 'SIE': sie_params, 'ANFW': nfw_params, 'PERT': pert_params}

def rms_extract(model_ver, model_path):
    """
    This is your original, detailed extraction logic, slightly adapted.
    """
    opt_result_file = os.path.join(model_path, f'{model_ver}_optresult.dat')
    
    with open(opt_result_file, 'r') as file:
        opt_result = file.readlines()

    last_optimize_index = next((i for i, line in reversed(list(enumerate(opt_result))) if 'optimize' in line), None)
    if last_optimize_index is None:
        raise ValueError("No 'optimize' line found in the file.")

    opt_result = opt_result[last_optimize_index + 1:]

    lens_params_dict = {}
    for line in opt_result:
        if line.startswith('lens'):
            parts = re.split(r'\s+', line.strip())
            lens_name = parts[1]
            params = [float(x) for x in parts[2:]]
            lens_params_dict[lens_name] = params[1:]

    source_params = []
    for line in opt_result:
        if line.startswith('point'):
            parts = re.split(r'\s+', line.strip())
            source_params.append([float(x) for x in parts[1:]])
    
    chi2_line = next((line for line in opt_result if 'chi^2' in line), None)
    if chi2_line is None:
        raise ValueError("No line with 'chi^2' found in the file.")
    chi2_value = float(chi2_line.split('=')[-1].strip().split()[0])
    
    hubble_val = None
    if h0:
        hubble_line = next((line for line in opt_result if 'hubble =' in line), None)
        if hubble_line:
            hubble_val = float(hubble_line.split('=')[-1].strip().split()[0])

    obs_point = pd.read_csv(obs_point_file, sep='\s+', header=None, skiprows=1, names=['x', 'y', 'mag', 'pos_err', 'mag_err', '1', '2', '3'])
    out_point_file = os.path.join(model_path, f'{model_ver}_point.dat')
    out_point = pd.read_csv(out_point_file, sep='\s+', header=None, skiprows=1, names=['x', 'y', 'mag', 'td', 'col5', 'col6', 'col7', 'col8'])
    
    td_vals = list(out_point['td']) if time_delay else None

    out_point['x_diff'] = abs(out_point['x'] - obs_point['x'])
    out_point['y_diff'] = abs(out_point['y'] - obs_point['y'])
    out_point['mag_diff'] = abs(abs(out_point['mag']) - abs(obs_point['mag']))
    out_point['pos_sq'] = np.sqrt((out_point['x_diff']**2 + out_point['y_diff']**2).astype(float))
    pos_rms = np.average(out_point['pos_sq'])
    mag_rms = np.average(np.sqrt((out_point['mag_diff']**2).astype(float)))
    
    return pos_rms, mag_rms, chi2_value, source_params, lens_params_dict, hubble_val, td_vals

def run_single_model(params, rank):
    m_val, n_val, o_val = params
    model_name = f'{model}_{m_val}_{n_val}_{o_val}'
    
    # --- MODIFIED: All temporary work is now done in the FAST LOCAL SCRATCH directory ---
    worker_temp_dir = os.path.join(SCRATCH_DIR, f'worker_temp_{rank}')
    os.makedirs(worker_temp_dir, exist_ok=True)
    temp_input_py_file = os.path.join(worker_temp_dir, f"temp_input_{uuid.uuid4()}.py")

    try:
        shutil.copy(input_py_file, temp_input_py_file)
        with open(temp_input_py_file, 'r') as f: lines = f.readlines()
        
        new_lines = []
        for line in lines:
            if line.strip().startswith('path ='):
                new_lines.append(f"path = '{os.path.join(worker_temp_dir, model_name)}'\n")
            elif line.strip().startswith('constraint_file ='):
                new_lines.append(f"constraint_file = '{constraint_file}'\n")
            elif prior_file and line.strip().startswith('prior_file ='):
                new_lines.append(f"prior_file = '{prior_file}'\n")
            else:
                new_lines.append(line)
        lines = new_lines

        lens_lines_indices = [i for i, line in enumerate(lines) if line.strip().startswith('glafic.set_lens(')]
        params_to_set = [(m_lens, m_param, m_val), (n_lens, n_param, n_val), (o_lens, o_param, o_val)]
        for lens_target, param_index, value in params_to_set:
            line_idx = lens_lines_indices[lens_target - 1]
            parts = [p.strip() for p in re.search(r'\((.*)\)', lines[line_idx]).group(1).split(',')]
            parts[2 + param_index - 1] = str(value)
            lines[line_idx] = re.sub(r'\(.*\)', f'({", ".join(parts)})', lines[line_idx])
        
        with open(temp_input_py_file, 'w') as f: f.writelines(lines)
        
        env = os.environ.copy()
        for key in ["MKL_NUM_THREADS", "OMP_NUM_THREADS"]: env[key] = "1"
        subprocess.run(['python3', temp_input_py_file], env=env, check=True, capture_output=True, text=True)
        
        pos_rms, mag_rms, chi2, source, lens_params, hubble_val, td_vals = rms_extract(model_name, worker_temp_dir)
        
        out_point_file = os.path.join(worker_temp_dir, f'{model_name}_point.dat')
        num_images = sum(1 for line in open(out_point_file) if line.strip()) - 1 if os.path.exists(out_point_file) else 0
        
        result_dict = {'m': m_val, 'n': n_val, 'o': o_val, 'num_images': num_images, 'pos_rms': pos_rms, 'mag_rms': mag_rms, 'chi2': chi2, 'source_x': source[0][1] if source else 0, 'source_y': source[0][2] if source else 0}
        if time_delay and td_vals is not None: result_dict['time_delays'] = ';'.join(map(str, td_vals))
        if h0 and hubble_val is not None: result_dict['h0'] = hubble_val

        for lens_name, params in lens_params.items():
            model_type = lens_name.upper()
            if model_type not in model_params: continue
            param_names = model_params[model_type]
            for i, param_val in enumerate(params):
                if i < len(param_names) and param_names[i] != 'NaN':
                    col_name = f"{model_type}_{param_names[i]}"
                    result_dict[col_name] = param_val

        return result_dict
    except Exception:
        return None
    finally:
        # Clean up the LOCAL temporary directory for this worker
        if os.path.exists(worker_temp_dir):
            shutil.rmtree(worker_temp_dir)

def main():
    rank = int(os.getenv('SLURM_PROCID', 0))
    size = int(os.getenv('SLURM_NTASKS', 1))
    
    from itertools import product
    all_tasks = [params for params in product(m, n, o)]
    tasks_for_this_worker = all_tasks[rank::size]
    
    output_csv_file = os.path.join(final_results_dir, f"results_rank_{rank}.csv")
    status_file = os.path.join(progress_status_dir, f"worker_{rank}.status")
    
    num_completed = 0
    total_for_this_worker = len(tasks_for_this_worker)
    
    for i, task_params in enumerate(tasks_for_this_worker):
        result = run_single_model(task_params, rank)
        if result:
            df = pd.DataFrame([result])
            header = not os.path.exists(output_csv_file)
            df.to_csv(output_csv_file, mode='a', header=header, index=False)
            num_completed += 1
            
        is_update_iteration = (i + 1) % PROGRESS_UPDATE_INTERVAL == 0
        is_last_iteration = (i + 1) == total_for_this_worker
        if is_update_iteration or is_last_iteration:
            with open(status_file, 'w') as f:
                f.write(str(i + 1))
            
    print(f"Worker {rank} finished. Total models attempted: {total_for_this_worker}. Successful: {num_completed}.")
    sys.stdout.flush()

if __name__ == "__main__":
    main()