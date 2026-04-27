import torch
import numpy as np
import time
from tqdm import tqdm

def train_loop(model, loader_train, loader_val, optimizer, criterion, num_epoch, device, scheduler=None):
    """
    Core training loop for N2V.
    
    Args:
        model: The PyTorch model (UNet).
        loader_train: DataLoader for training data.
        loader_val: DataLoader for validation data.
        optimizer: The optimizer.
        criterion: The loss function.
        num_epoch: Number of epochs to train.
        device: 'cuda' or 'cpu'.
        scheduler: Learning rate scheduler (optional).
        
    Returns:
        model: The trained model.
        history: Dictionary containing training history (loss, etc.).
    """
    model = model.to(device)
    
    history = {
        'train_loss': [],
        'val_loss': []
    }
    
    print(f"Starting training on {device} for {num_epoch} epochs...")
    
    for epoch in range(1, num_epoch + 1):
        start_time = time.time()
        
        # --- Training Phase ---
        model.train()
        epoch_loss = 0.0
        num_batches = 0
        
        with tqdm(total=len(loader_train), desc=f'Epoch {epoch}/{num_epoch}', unit='batch') as pbar:
            for batch_idx, data in enumerate(loader_train):
                # Data is a dict from MemoryDataset
                # input: masked image
                # label: original noisy image
                # mask: 0.0 for masked pixels, 1.0 for others
                
                input_img = data['input'].to(device)
                label_img = data['label'].to(device)
                mask = data['mask'].to(device)
                
                optimizer.zero_grad()
                
                output = model(input_img)
                
                # N2V Loss: Calculate loss only on masked pixels (where mask == 0.0)
                # The formula (1 - mask) ensures we only keep values where mask is 0.0
                loss = criterion(output * (1 - mask), label_img * (1 - mask))
                
                loss.backward()
                optimizer.step()
                
                epoch_loss += loss.item()
                num_batches += 1
                
                pbar.set_postfix({'loss': loss.item()})
                pbar.update(1)
        
        avg_train_loss = epoch_loss / num_batches if num_batches > 0 else 0.0
        history['train_loss'].append(avg_train_loss)
        
        # --- Validation Phase ---
        avg_val_loss = 0.0
        if loader_val is not None and len(loader_val) > 0:
            model.eval()
            val_loss = 0.0
            val_batches = 0
            
            with torch.no_grad():
                for data in loader_val:
                    input_img = data['input'].to(device)
                    label_img = data['label'].to(device)
                    mask = data['mask'].to(device)
                    
                    output = model(input_img)
                    loss = criterion(output * (1 - mask), label_img * (1 - mask))
                    
                    val_loss += loss.item()
                    val_batches += 1
            
            avg_val_loss = val_loss / val_batches if val_batches > 0 else 0.0
            history['val_loss'].append(avg_val_loss)
            
            if scheduler:
                if isinstance(scheduler, torch.optim.lr_scheduler.ReduceLROnPlateau):
                    scheduler.step(avg_val_loss)
                else:
                    scheduler.step()
                    
            print(f"Epoch {epoch} | Train Loss: {avg_train_loss:.6f} | Val Loss: {avg_val_loss:.6f} | Time: {time.time() - start_time:.2f}s")
        else:
            print(f"Epoch {epoch} | Train Loss: {avg_train_loss:.6f} | Time: {time.time() - start_time:.2f}s")
            
    print("Training finished.")
    return model, history
