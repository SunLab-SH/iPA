# Use PyTorch official image with CUDA support (Python 3.9)
# If GPU is not available, PyTorch will automatically fallback to CPU
FROM pytorch/pytorch:1.13.1-cuda11.6-cudnn8-runtime

# Set working directory
WORKDIR /app

# Install system dependencies required by some Python packages (e.g., OpenCV, imageio)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgl1-mesa-glx \
    libglib2.0-0 \
    libsm6 \
    libxext6 \
    libxrender-dev \
    git-lfs \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements first to leverage Docker cache
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy the iPA source code
COPY ipa/ ./ipa/
COPY setup.py .
COPY README.md .
COPY MANIFEST.in .

# Install iPA in development mode
RUN pip install -e .

# Initialize Git LFS (to pull model weights if needed during build, 
# but usually models are pulled at runtime or via volume)
RUN git lfs install

# Create a non-root user for security (optional but recommended)
RUN useradd -m -u 1000 ipa_user && chown -R ipa_user:ipa_user /app
USER ipa_user

# Default command
CMD ["/bin/bash"]
