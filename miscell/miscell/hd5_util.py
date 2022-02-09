import h5py
import numpy as np

def save_d_result_to_hd5(h5_name, columns, d_result):
    h5_handle = h5py.File(h5_name, 'w') 
    for key in columns:
        h5_handle.create_dataset(key, data=np.array(d_result[key]))
    h5_handle.close()
    print(f'make {h5_name}')

def read_d_result_from_hd5(h5_name, columns):
    h5_handle = h5py.File(h5_name, 'r')
    d_result = {key: np.array(h5_handle.get(f'{key}')) for key in columns}
    h5_handle.close()
    print(f'Read {h5_name}')
    return d_result
