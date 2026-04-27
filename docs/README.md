# iPA Documentation

This directory contains the Sphinx documentation for the iPA (Integrated Processing and Analysis) toolkit.

## Building the Documentation

### Prerequisites

Make sure you have Python 3.7+ installed, then install the required packages:

```bash
pip install -r requirements.txt
```

### Building HTML Documentation

On Windows:

```cmd
make.bat html
```

On Linux/macOS:

```bash
make html
```

The generated HTML documentation will be available in `_build/html/index.html`.

### Building PDF Documentation (Optional)

To build PDF documentation, you'll need LaTeX installed:

```bash
make latexpdf
```

## Documentation Structure

- `index.rst` - Main documentation entry point
- `conf.py` - Sphinx configuration file
- `rest_files/` - Individual module documentation
  - `processing.rst` - Processing module documentation
  - `et_analysising.rst` - CryoET analysis module documentation
  - `partitioning.rst` - Spatial partitioning module documentation
  - `common.rst` - Common utilities documentation
  - `examples.rst` - Examples and tutorials
- `_static/` - Static files (CSS, images, etc.)
- `_templates/` - Custom Sphinx templates

## Customization

- Edit `conf.py` to modify Sphinx settings
- Modify `_static/custom.css` for visual customization
- Add images to `rest_files/resources/` directory
- Update individual `.rst` files to modify content

## API Documentation

The documentation uses Sphinx autodoc to automatically generate API documentation from docstrings in the source code. Make sure your Python modules are properly documented with docstrings following the NumPy/Google docstring format.

## Contributing

When adding new modules or functions to iPA:

1. Ensure proper docstrings are included
2. Add corresponding `.rst` files if needed
3. Update the main `index.rst` toctree
4. Rebuild documentation to verify changes

For questions about documentation, contact: <liad@shanghaitech.edu.cn>
