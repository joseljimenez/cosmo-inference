import logging
logger = logging.getLogger(__name__)
import bin_range
from pathlib import Path
from cosmology.cosmology import Cosmology
from astropy.io import fits
import numpy as np

a = bin_range.zmin
b = bin_range.zmax

try:
    SCRIPT_DIR = Path(__file__).resolve().parent
except:
    SCRIPT_DIR = Path.cwd() 
BASE_DIR = SCRIPT_DIR.parent

#-----------------------------------------------------------------------
# ---------------------------- UNION3  ---------------------------------
#-----------------------------------------------------------------------
# Import the data 

# file_Union3   = fits.open(BASE_DIR+"/data/union3_release/mu_mat_union3_cosmo=2_mu.fits")
file_Union3 = fits.open(
    BASE_DIR / "data" / "union3_release" / "mu_mat_union3_cosmo=2_mu.fits"
)
data_Union3   = file_Union3[0].data
zcmb_Union3   = data_Union3[0,1 :]                           # Redshift corrected for CMB
# zhel_Union3   = data_Union3[0, 1:]                         # Heliocentric redshift 
mb_Union3     = data_Union3[1: ,0]                           # data of the apparent magnitude mb
cov_matUnion3 = np.linalg.inv(data_Union3[1:,1:])           # Covariance matrix

mask_Union3 = (zcmb_Union3<=b) & (zcmb_Union3>a)                   # mask to separate data

# Bin the data using the mask

zcmb_bin_Union3   = zcmb_Union3[mask_Union3]
mb_bin_Union3     = mb_Union3[mask_Union3]
covbin_mat_Union3 = cov_matUnion3[mask_Union3, :][:, mask_Union3]     # Reconstruct the covariance matrix according to the binned data
inv_cov_Union3    = np.linalg.inv(covbin_mat_Union3)                  # Inverse of the covariance matrix

# Theoretical apparent magnitude 
def mb_theo_Union3(Ob0,H0,As,alpha,M):
    cosmo = Cosmology(Ob0, H0, As, alpha)
    return np.array([cosmo.mb(zi, M) for zi in zcmb_bin_Union3])   

# NOTE: In the Union3 and DES samples, the fitted magnitude parameter represents 
# an offset ΔM relative to the absolute magnitude calibration (MB) rather than MB itself. 
# Therefore, M here should be interpreted as ΔM, not as the absolute magnitude.

# Construct the chi^2 for SNe - Union3 
def chi2_Union3(Ob0,H0,As,alpha,M):
    delta = mb_bin_Union3 - mb_theo_Union3(Ob0, H0, As, alpha, M)
    chi2  = np.dot(delta, np.dot(inv_cov_Union3, delta))
    return chi2

#-----------------------------------------------------------------------
# ---------------------------- DESI-BAO DR2-----------------------------
#-----------------------------------------------------------------------
# DESI BAO data
arr = np.genfromtxt(
    BASE_DIR / "data" / "BAO_data" / "desi_bao_dr2" /
    "desi_gaussian_bao_ALL_GCcomb_mean.txt",
    skip_header=1, dtype=None, encoding=None
)

z_desi    = np.array([d[0] for d in arr], dtype=float)
data_desi = np.array([d[1] for d in arr], dtype=float)
qty_desi  = np.array([d[2] for d in arr])

# Covariance matrix
covmat_desi = np.loadtxt(
    BASE_DIR / "data" / "BAO_data" / "desi_bao_dr2" /
    "desi_gaussian_bao_ALL_GCcomb_cov.txt"
)

# Redshift mask
mask_DESI = (z_desi <= b) & (z_desi > a)

z_desi_sel    = z_desi[mask_DESI]
data_desi_sel = data_desi[mask_DESI]
qty_desi_sel  = qty_desi[mask_DESI]
cov_desi_sel  = covmat_desi[mask_DESI, :][:, mask_DESI]
inv_cov_desi  = np.linalg.inv(cov_desi_sel)

def DM(z, Ob0, H0, As, alpha):
    cosmo = Cosmology(Ob0, H0, As, alpha)
    return np.array([cosmo.DM(zi) for zi in z])

def DH(z, Ob0, H0, As, alpha):
    cosmo = Cosmology(Ob0, H0, As, alpha)
    return np.array([cosmo.DH(zi) for zi in z])

def chi2_DESI(Ob0, H0, As, alpha, rd):
    DMv = DM(z_desi_sel, Ob0, H0, As, alpha)
    DHv = DH(z_desi_sel, Ob0, H0, As, alpha)

    model = np.empty_like(data_desi_sel)

    for i in range(len(z_desi_sel)):
        if qty_desi_sel[i] == "DM_over_rs":
            model[i] = DMv[i] / rd
        elif qty_desi_sel[i] == "DH_over_rs":
            model[i] = DHv[i] / rd
        elif qty_desi_sel[i] == "DV_over_rs":
            model[i] = (z_desi_sel[i] * DHv[i] * DMv[i]**2)**(1.0 / 3.0) / rd
        else:
            raise ValueError(f"Cantidad DESI no reconocida: {qty_desi_sel[i]}")

    r = data_desi_sel - model
    return float(r @ inv_cov_desi @ r)


#-----------------------------------------------------------------------
# ----------------------------- BBN ------------------------------------
#-----------------------------------------------------------------------
BBN_mean = 2.166
BBN_sigma = np.sqrt(0.015**2 + 0.011**2)

def chi2_BBN(Ob0, H0):
    h = H0 / 100.0
    Obh2 = 100.0 * Ob0 * h**2
    chi2 = (Obh2 - BBN_mean)**2 / BBN_sigma**2
    logger.debug(
        f"BBN prior | Ob0={Ob0:.5e}, H0={H0:.3f}, "
        f"h={h:.5f}, 100*Ob0*h^2={Obh2:.5f}, "
        f"chi2_BBN={chi2:.5f}"
    )
    return chi2


def build_total_chi2():
    param_names = ["Ob0", "H0", "As", "alpha", "M", "rd"]
    def chi2_total(theta):
        Ob0, H0, As, alpha, M, rd = theta
        chi2 = 0.0
        # SNe Ia (Union3)
        chi2 += chi2_Union3(Ob0, H0, As, alpha, M)
        # BBN prior
        chi2 += chi2_BBN(Ob0, H0)
        # DESI BAO
        chi2 += chi2_DESI(Ob0, H0, As, alpha, rd)

        return float(chi2)

    return chi2_total, param_names
    
