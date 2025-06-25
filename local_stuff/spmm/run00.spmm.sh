
PREV_DIR=$(readlink -f .)
BUILD_DIR="${PREV_DIR}/../../cmake-build-debug"
#BUILD_DIR="${PREV_DIR}/../../build"

DATA_DIR="/Users/peng599/Library/CloudStorage/OneDrive-PNNL/Documents/Datasets"

OUTPUT_DIR="output.$(date +%FT%T)"
mkdir -p ${OUTPUT_DIR}

cd ${OUTPUT_DIR}

dram_spmm_app="${BUILD_DIR}/run_dram_spmm_csr"
#rapid_gemm_app="${BUILD_DIR}/run_rapid_spmm_csr"

# Run
${dram_spmm_app} ${DATA_DIR}
#${rapid_gemm_app} ${DATA_DIR}

