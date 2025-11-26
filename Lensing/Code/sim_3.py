#!/usr/bin/env python
import os
import traceback
import numpy as np
import pandas as pd
import re
import sys
import shutil
import uuid
from itertools import product, islice

# IMPORTANT: Ensure the 'glafic' library is in your Python path
import glafic

# --- File Paths & Parameters ---
model_output_base = '/home/rommulus/Projects/itng_lensing/Simulations/Output/Sim_3'
base_results_path = '/home/rommulus/Projects/itng_lensing/Simulations/Input/System_2'
obs_point_file = os.path.join(base_results_path, 'pos+flux_point.dat')
SCRATCH_DIR = os.getenv('SCRATCH_DIR', '/tmp')

# --- Simulation Parameters ---
sim_name = 'Sim 23'
model = 'NFW'
m = np.linspace(0.001, 0.1, 1000)
n = np.linspace(0, 360, 1000)
o = np.array([0.0])

m_lens, m_param = 2, 5
n_lens, n_param = 2, 6
o_lens, o_param = 2, 8

constraint_file = os.path.join(base_results_path, 'pos_point.dat')
prior_file = None
time_delay = False
h0 = False
critical_curve = False

# --- Directory setup ---
final_results_dir = os.path.join(model_output_base, 'individual_results')
progress_status_dir = os.path.join(model_output_base, 'progress_status')
os.makedirs(final_results_dir, exist_ok=True)
os.makedirs(progress_status_dir, exist_ok=True)

# --- Performance Tuning ---
PROGRESS_UPDATE_INTERVAL = 100
CSV_BATCH_SIZE = 100

# --- Parameter definitions (used for CSV headers) ---
pow_params = ['$z_{s,fid}$', 'x', 'y', 'e', '$θ_{e}$', '$r_{Ein}$', '$\gamma$ (PWI)']
sie_params = ['$\sigma$', 'x', 'y', 'e', '$θ_{e}$', '$r_{core}$', 'NaN']
nfw_params = ['M', 'x', 'y', 'e', '$θ_{e}$', 'c or $r_{s}$', 'NaN']
pert_params = ['$z_{s,fid}$', 'x', 'y', '$\gamma$', '$θ_{\gamma}$', 'NaN', '$\kappa$']
model_params = {'POW': pow_params, 'SIE': sie_params, 'ANFW': nfw_params, 'PERT': pert_params}


# --- Helper Functions ---
def calculate_distance(row1, row2):
    return np.sqrt((row1['x'] - row2['x'])**2 + (row1['y'] - row2['y'])**2)

def assign_image_names(df, start_index):
    names = ['A', 'B', 'C', 'D']
    current_name_index = 0
    current_index = start_index
    assigned_indices = {current_index}
    
    while len(assigned_indices) < len(df):
        current_x, current_y = df.at[current_index, 'x'], df.at[current_index, 'y']
        next_index = None
        max_angle = -float('inf')
        
        for i in range(len(df)):
            if i in assigned_indices:
                continue
            dx = df.at[i, 'x'] - current_x
            dy = df.at[i, 'y'] - current_y
            angle = np.arctan2(dy, dx)
            if angle > max_angle:
                max_angle = angle
                next_index = i
        
        if next_index is not None:
            current_name_index += 1
            df.at[next_index, 'Img'] = names[current_name_index]
            assigned_indices.add(next_index)
            current_index = next_index

obs_point = pd.read_csv(obs_point_file, sep='\s+', header=None, skiprows=1, names=['x', 'y', 'mag', 'pos_err', 'mag_err', 'td', 'td_err', 'parity'])
brightest_index = obs_point['mag'].idxmax()
obs_point.at[brightest_index, 'Img'] = 'A'
assign_image_names(obs_point, brightest_index)
n_img = len(obs_point)

def rms_extract(model_ver, model_path, obs_point=obs_point):
    opt_result_file = os.path.join(model_path, f'{model_ver}_optresult.dat')
    with open(opt_result_file, 'r') as file: opt_result = file.readlines()
    last_optimize_index = next((i for i, line in reversed(list(enumerate(opt_result))) if 'optimize' in line), None)
    if last_optimize_index is None: raise ValueError("No 'optimize' line found.")
    opt_result = opt_result[last_optimize_index + 1:]
    lens_params_dict = {}
    for line in opt_result:
        if line.startswith('lens'):
            parts = re.split(r'\s+', line.strip())
            lens_params_dict[parts[1]] = [float(x) for x in parts[3:]]
    source_params = [[float(x) for x in re.split(r'\s+', line.strip())[1:]] for line in opt_result if line.startswith('point')]
    chi2_line = next((line for line in opt_result if 'chi^2' in line), None)
    if chi2_line is None: raise ValueError("No 'chi^2' line found.")
    chi2_value = float(chi2_line.split('=')[-1].strip().split()[0])
    hubble_val = None
    if h0:
        hubble_line = next((line for line in opt_result if 'hubble =' in line), None)
        if hubble_line: hubble_val = float(hubble_line.split('=')[-1].strip().split()[0])
    out_point_file = os.path.join(model_path, f'{model_ver}_point.dat')
    out_point = pd.read_csv(out_point_file, sep='\s+', header=None, skiprows=1, names=['x', 'y', 'mag', 'td', 'col5', 'col6', 'col7', 'col8'])
    num_pred_images = len(out_point)

    if num_pred_images > n_img:
        distance_matrix = np.zeros((n_img, num_pred_images))
        for i in range(n_img):
            for j in range(num_pred_images):
                distance_matrix[i, j] = calculate_distance(obs_point.iloc[i], out_point.iloc[j])

        matched_indices = np.unravel_index(np.argsort(distance_matrix, axis=None), distance_matrix.shape)
        matched_pred_indices = set()
        matched_obs_indices = set()
        matches = []
        for obs_idx, pred_idx in zip(*matched_indices):
            if obs_idx not in matched_obs_indices and pred_idx not in matched_pred_indices:
                matches.append((obs_idx, pred_idx))
                matched_obs_indices.add(obs_idx)
                matched_pred_indices.add(pred_idx)
            if len(matches) == n_img:
                break

        matches.sort(key=lambda x: x[0])
        matched_pred_indices = [pred_idx for _, pred_idx in matches]
        out_point = out_point.iloc[matched_pred_indices].reset_index(drop=True)
        out_point['Img'] = obs_point['Img'].values
    
    if num_pred_images < n_img:
        distance_matrix = np.zeros((n_img, num_pred_images))
        for i in range(n_img):
            for j in range(num_pred_images):
                distance_matrix[i, j] = calculate_distance(obs_point.iloc[i], out_point.iloc[j])

        matched_indices = np.unravel_index(np.argsort(distance_matrix, axis=None), distance_matrix.shape)
        matched_pred_indices = set()
        matched_obs_indices = set()
        matches = []
        for obs_idx, pred_idx in zip(*matched_indices):
            if obs_idx not in matched_obs_indices and pred_idx not in matched_pred_indices:
                matches.append((obs_idx, pred_idx))
                matched_obs_indices.add(obs_idx)
                matched_pred_indices.add(pred_idx)
            if len(matches) == num_pred_images:
                break

        matches.sort(key=lambda x: x[0])
        matched_pred_indices = [pred_idx for _, pred_idx in matches]
        out_point = out_point.iloc[matched_pred_indices].reset_index(drop=True)
        out_point['Img'] = [obs_point.at[obs_idx, 'Img'] for obs_idx, _ in matches]

    if num_pred_images == n_img:
        distance_matrix = np.zeros((n_img, n_img))
        for i in range(n_img):
            for j in range(n_img):
                distance_matrix[i, j] = calculate_distance(obs_point.iloc[i], out_point.iloc[j])
        
        matched_indices = np.unravel_index(np.argsort(distance_matrix, axis=None), distance_matrix.shape)
        matched_pred_indices = set()
        matched_obs_indices = set()
        matches = []
        for obs_idx, pred_idx in zip(*matched_indices):
            if obs_idx not in matched_obs_indices and pred_idx not in matched_pred_indices:
                matches.append((obs_idx, pred_idx))
                matched_obs_indices.add(obs_idx)
                matched_pred_indices.add(pred_idx)
            if len(matches) == n_img:
                break
        matches.sort(key=lambda x: x[0])
        matched_pred_indices = [pred_idx for _, pred_idx in matches]
        out_point = out_point.iloc[matched_pred_indices].reset_index(drop=True)
        out_point['Img'] = [obs_point.at[obs_idx, 'Img'] for obs_idx, _ in matches]
    
    obs_point = obs_point.sort_values(by='Img').reset_index(drop=True)
    out_point = out_point.sort_values(by='Img').reset_index(drop=True)

    image_rms = []
    pos_rms = np.nan

    num_matched_images = len(out_point)

    if num_matched_images > 0:
        for i in range(num_matched_images):
            obs_row = obs_point[obs_point['Img'] == out_point.at[i, 'Img']]
            if not obs_row.empty:
                dist = np.sqrt((obs_row.iloc[0]['x'] - out_point.at[i, 'x'])**2 + 
                               (obs_row.iloc[0]['y'] - out_point.at[i, 'y'])**2)
                image_rms.append(dist)

        if image_rms:
            pos_rms = np.sqrt(np.sum(np.array(image_rms)**2) / num_matched_images)

    image_rms = np.array(image_rms)

    flux_rms = []
    mag_rms = np.nan
    for i in range(len(out_point)):
        diff = abs(abs(out_point.at[i, 'mag']) - abs(obs_point.at[i, 'mag']))
        flux_rms.append(diff)
        mag_rms = np.sqrt(np.sum(np.array(flux_rms)**2) / len(out_point))
    flux_rms = np.array(flux_rms)

    # Percentage Errors in Predicted Magnification/Flux
    percentage_errors = []
    for i in range(len(out_point)):
        perc_error = abs(abs((out_point.at[i, 'mag']) - abs(obs_point.at[i, 'mag']))) / abs(obs_point.at[i, 'mag']) * 100
        percentage_errors.append(perc_error)
    avg_percentage_error = np.mean(percentage_errors) if len(percentage_errors) > 0 else 0
    percentage_errors = np.array(percentage_errors)

    percentage_errors_td = []
    avg_percentage_error_td = None
    td_rms = None
    if time_delay:
            for i in range(len(out_point)):
                diff = out_point.at[i, 'td'] - obs_point.at[i, 'td']
                out_point.at[i, 'td_rms'] = diff
            percentage_errors_td = [abs(out_point.at[i, 'td_rms'] / obs_point.at[i, 'td']) * 100 for i in range(len(out_point)) if obs_point.at[i, 'td'] != 0]
            for i in range(len(out_point)):
                if obs_point.at[i, 'td'] == 0:
                    percentage_errors_td.insert(i, 0)
            avg_percentage_error_td = np.mean(percentage_errors_td) if len(percentage_errors_td) > 0 else 0
            percentage_errors_td = np.array(percentage_errors_td)
            td_rms = np.sqrt(np.sum(out_point['td_rms']**2) / len(out_point))
    else:
        td_rms = None
        percentage_errors_td = None
        avg_percentage_error_td = None
    
    td_vals = np.array(out_point['td']) if time_delay else None
    return pos_rms, image_rms, mag_rms, flux_rms, percentage_errors, avg_percentage_error, chi2_value, source_params, lens_params_dict, hubble_val, td_vals, td_rms, percentage_errors_td, avg_percentage_error_td, out_point


def run_glafic_calculation(params, model_name, worker_temp_dir):
    m_val, n_val, o_val = params
    output_path = os.path.join(worker_temp_dir, model_name)
    
    base_lens_params = [0.261343256161012, 1e11, 0.0, 0.0, 0.2, 0.0, 80, 0.0]
    base_shear_params = [0.261343256161012, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    current_lens_params = list(base_lens_params)
    current_shear_params = list(base_shear_params)
    current_shear_params[n_param - 1] = n_val
    current_shear_params[o_param - 1] = o_val
    current_shear_params[m_param - 1] = m_val

    glafic.init(0.3089901684739047, 0.6910098315260953, -1.0, 0.6736, output_path, -3.0, -3.0, 3.0, 3.0, 0.01, 0.01, 1, verb=0)
    glafic.set_secondary('chi2_splane 1', verb=0)
    glafic.set_secondary('chi2_checknimg 0', verb=0)
    glafic.set_secondary('chi2_restart   -1', verb=0)
    glafic.set_secondary('chi2_usemag    1', verb=0)
    glafic.set_secondary('hvary          0', verb=0)
    glafic.set_secondary('ran_seed -122000', verb=0)
    glafic.startup_setnum(2, 0, 1)
    glafic.set_lens(1, 'anfw', *current_lens_params)
    glafic.set_lens(2, 'pert', *current_shear_params)
    glafic.set_point(1, 1.0, 0.0, 0.0)
    glafic.setopt_lens(1, 0, 1, 1, 1, 1, 1, 1, 0)
    glafic.setopt_lens(2, 0, 0, 0, 0, 1, 1, 0, 0)
    glafic.setopt_point(1, 0, 1, 1)
    glafic.model_init(verb=0)
    glafic.readobs_point(constraint_file)
    if prior_file and os.path.exists(prior_file):
        glafic.parprior(prior_file)
    glafic.optimize()
    glafic.findimg()
    
    if critical_curve:
        glafic.writecrit(1.0)
        
    glafic.quit()


def run_single_model(params, worker_temp_dir, obs_point_df):
    m_val, n_val, o_val = round(params[0], 5), round(params[1], 5), round(params[2], 5)
    model_name = f'{model}_{m_val}_{n_val}_{o_val}'

    try:
        # Redirect both Python-level and OS-level stdout/stderr to devnull
        devnull = open(os.devnull, 'w')
        _saved_stdout_fd = os.dup(1)
        _saved_stderr_fd = os.dup(2)
        os.dup2(devnull.fileno(), 1)
        os.dup2(devnull.fileno(), 2)
        _saved_stdout = sys.stdout
        _saved_stderr = sys.stderr
        sys.stdout = devnull
        sys.stderr = devnull
        run_glafic_calculation((m_val, n_val, o_val), model_name, worker_temp_dir)
        
        pos_rms, image_rms, mag_rms, flux_rms, percentage_errors, avg_percentage_error, chi2, source, lens_params, hubble_val, td_vals, td_rms, percentage_errors_td, avg_percentage_error_td, out_point = rms_extract(model_name, worker_temp_dir)
        
        out_point_file = os.path.join(worker_temp_dir, f'{model_name}_point.dat')
        num_images = sum(1 for line in open(out_point_file) if line.strip()) - 1 if os.path.exists(out_point_file) else 0
        
        result_dict = {'m': m_val, 'n': n_val, 'o': o_val, 'num_images': num_images, 'pos_rms': pos_rms, 'mag_rms': mag_rms, 'avg_mag_per': avg_percentage_error, 'chi2': chi2, 'source_x': source[0][1] if source else 0, 'source_y': source[0][2] if source else 0}
        for i, row in obs_point.iterrows():
            image_name = row['Img']
            result_dict[f'x_{image_name}'] = out_point.at[i, 'x'] if i < len(out_point) else None
            result_dict[f'y_{image_name}'] = out_point.at[i, 'y'] if i < len(out_point) else None
            result_dict[f'pos_rms_{image_name}'] = image_rms[i] if i < len(image_rms) else None
            result_dict[f'mag_rms_{image_name}'] = flux_rms[i] if i < len(flux_rms) else None
            result_dict[f'mag_per_{image_name}'] = percentage_errors[i] if i < len(percentage_errors) else None
        if time_delay and td_vals is not None:
            for i, row in obs_point.iterrows():
                image_name = row['Img']
                result_dict[f'td_{image_name}'] = td_vals[i] if i < len(td_vals) else None
                result_dict[f'td_rms_{image_name}'] = out_point.at[i, 'td_rms'] if i < len(out_point) else None
                result_dict[f'td_per_{image_name}'] = percentage_errors_td[i] if i < len(percentage_errors_td) else None

        if time_delay and avg_percentage_error_td is not None: result_dict['avg_td_per'] = avg_percentage_error_td
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
        print(f"!!! WORKER FAILED on model {model_name}.")
        traceback.print_exc()
        sys.stdout.flush()
        return None
    finally:
        for suffix in ['_optresult.dat', '_point.dat', '_crit.dat']:
             f = os.path.join(worker_temp_dir, f"{model_name}{suffix}")
             if os.path.exists(f):
                 try: os.remove(f)
                 except OSError: pass

def write_batch_to_csv(batch, csv_file):
    if not batch: return
    df = pd.DataFrame(batch)
    header = not os.path.exists(csv_file)
    df.to_csv(csv_file, mode='a', header=header, index=False)

def main():
    rank = int(os.getenv('SLURM_PROCID', 0))
    size = int(os.getenv('SLURM_NTASKS', 1))
    
    all_tasks_generator = product(m, n, o)
    tasks_for_this_worker = islice(all_tasks_generator, rank, None, size)
    
    total_combinations = len(m) * len(n) * len(o)
    total_for_this_worker = (total_combinations - rank + size - 1) // size
    
    output_csv_file = os.path.join(final_results_dir, f"results_rank_{rank}.csv")
    status_file = os.path.join(progress_status_dir, f"worker_{rank}.status")
    
    worker_temp_dir = os.path.join(SCRATCH_DIR, f'worker_temp_{rank}_{uuid.uuid4()}')
    os.makedirs(worker_temp_dir, exist_ok=True)
    
    results_batch = []
    num_completed = 0
    
    try:
        for i, task_params in enumerate(tasks_for_this_worker):
            result = run_single_model(task_params, worker_temp_dir, obs_point_df = obs_point)
            
            if result:
                results_batch.append(result)
                num_completed += 1
            
            if len(results_batch) >= CSV_BATCH_SIZE or ((i + 1) == total_for_this_worker and results_batch):
                write_batch_to_csv(results_batch, output_csv_file)
                results_batch = [] 

            if (i + 1) % PROGRESS_UPDATE_INTERVAL == 0 or (i + 1) == total_for_this_worker:
                with open(status_file, 'w') as f:
                    f.write(str(i + 1))
    finally:
        if os.path.exists(worker_temp_dir):
            shutil.rmtree(worker_temp_dir)

    print(f"Worker {rank} finished. Total models attempted: {total_for_this_worker}. Successful: {num_completed}.")
    sys.stdout.flush()

if __name__ == "__main__":
    main()