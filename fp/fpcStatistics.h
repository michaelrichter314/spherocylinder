#ifndef fpcstatistics_h
#define fpcstatistics_h

#ifndef FPC_BKS
    #error Need to define FPC_BKS (number of blocks)
#endif

#ifndef FPC_TPB
    #error Need to define FPC_TPB (threads per block)
#endif



namespace fp
{

    template <typename type>
	type max(type* data, int N)
	{
            __global__ void kernelMax(type* kData, type* kOut)
            {
                __shared__ type kTemp[FPC_TPB];
                
                int start = threadIdx.x + blockIdx.x*BlockDim.x;
                int end   = 
                
            }
            
            // first, 
            if (N > 0)
            {
                    
                    type max = data[0];
                    for (int i = 1; i < N; i++)
                            if (data[i] > max) {max = data[i];}
                    
                    return max;
            }
	}

}


#endif