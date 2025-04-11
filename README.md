# 📦 SVDGP : Singular Value Decomposition and Gaussian process modelling

**SVDGP** package provides tools for the decomposition, modeling, and prediction of spatio-temporal functions represented as:  
  
  $$f : (x, t) \in \mathbb{R}^{d} \times \mathbb{R} \mapsto f(x, t) \in \mathbb{R},$$
  
where the function $f(x, t)$ is allways observed for a spatial locations $x \in \mathbb{R}^{d}$ at discrete time points $t_{1}, \ldots, t_{N_t}$.

## 📥 Installation

You can install the latest version of the package manually or directly from GitHub.

### Option 1: Install from GitHub (Recommended)

Make sure you have the `devtools` package installed, then use:

```r
install.packages("devtools")
devtools::install_github("TheseAdama/SVDGP")
```
### Option 2: Manual Installation (Download & Install ZIP)

1. **Download the ZIP or TAR.GZ file**
   Download the latest version of the package in ZIP or TAR.GZ format.

   - For **Windows**: `SVDGP_x.y.z.zip`
   - For **Linux/macOS**: `SVDGP_x.y.z.tar.gz`

2. **Install the package manually in R**

   Open your R session and run one of the following commands, replacing the file path with where you downloaded the archive:

   - **On Windows**:
     ```r
     install.packages("path/to/SVDGP_x.y.z.zip", repos = NULL, type = "win.binary")
     ```

   - **On Linux/macOS**:
     ```r
     install.packages("path/to/SVDGP_x.y.z.tar.gz", repos = NULL, type = "source")
     ```

---

After installation, load the package:

```r
library(SVDGP)
```
