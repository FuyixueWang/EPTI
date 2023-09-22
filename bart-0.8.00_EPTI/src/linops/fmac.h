
#include <complex.h>

struct linop_s;
extern const struct linop_s* linop_fmac_create(unsigned int N, const long dims[N], 
		unsigned int oflags, unsigned int iflags, unsigned int flags, const complex float* tensor);

const struct linop_s* linop_PartitionDim_create(unsigned int N, const long dims[N], 
		const long Dim1, const long Dim2, const long K);

void set_fmac_tensor(const struct linop_s* InLinop, complex float* tensor);

struct fmac_data;
void setfmacdataToNewTensor(struct fmac_data * data, const complex float* tensor);
struct fmac_data * getpdata();

extern const struct linop_s* linop_fmacOnCPU_create(unsigned int N, const long dims[N], 
		unsigned int oflags, unsigned int iflags, unsigned int flags, const complex float* tensor);

extern void linop_fmac_set_tensor(const struct linop_s* lop, int N, const long tdims[N], const complex float* tensor);
