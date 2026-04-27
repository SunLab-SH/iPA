
import numpy as np
import mrcfile
import os

def create_dummy_data(shape=(32, 256, 256)):
    print("Creating dummy data...")
    clean = np.zeros(shape, dtype=np.float32)
    # Add some structure
    for i in range(shape[0]):
        clean[i, 50:200, 50:200] = 0.8 * np.sin(i/5)
    
    noise = np.random.normal(0, 0.1, shape).astype(np.float32)
    noisy = clean + noise
    return np.clip(noisy, 0, 1)

output_path = os.path.join(os.path.dirname(__file__), 'dummy_input.mrc')
data = create_dummy_data()
with mrcfile.new(output_path, overwrite=True) as mrc:
    mrc.set_data(data)
print(f"Saved {output_path}")
