#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include<time.h>
#include<stdbool.h>

/*
Date : 28th July 2022
Description: This code performs polynomial multiplication using
Iterative Radix2 NTT and Inverse NTT for two arrays of power of 2.
This is general purpose as it generates the mod and root for every iteration of arrays.

Check the website https://www.nayuki.io/page/number-theoretic-transform-integer-dft
for a detailed explanation of the code*/

long maxe(long x, long y)
{
    if(x>=y)
        return x;
    else
        return y;
}
long power(long x, long y, long mod)
{
    long z=1;
    for(long i=0;i<y;i++)
        z=((z%mod)*(x%mod))%mod;
    return z;
}
bool is_prime(long n)
{

    long flag=1;
    for (long i = 2; i <= sqrt(n); i++) {
 
        if (n % i == 0) {
            flag = 0;
            break;
        }
    }
 
    if (n <= 1)
        flag = 0;
 
    if (flag == 1)
        return true;
    else
        return false;
        
}


long * unique_prime_factors(long n)
{
    long* prime = (long*)malloc(n * sizeof(long)),sz=0;
    long i=2,end=sqrt(n);

    if(n<1)
        return NULL;
    while(i<=end)
    {
        if(n%i==0)
        {
            n=n/i;
            prime[sz]=i;
            sz++;
            while(n%i == 0)
                n=n/i;
            end=sqrt(n);
        }
        i+=1;
    }
    if(n>1)
    {
        prime[sz]=n;
        sz++;
    }
    return prime; 
}

long reciprocal(long a, long m)
{
    long m0 = m;
    long y = 0, x = 1;
 
    if (m == 1)
        return 0;
 
    while (a > 1) {
       
        long q = a / m;
        long t = m;
 
        m = a % m, a = t;
        t = y;
 
        
        y = x - q * y;
        x = t;
    }
 
    
    if (x < 0)
        x += m0;
 
    return x;
}
bool is_primitive(long val, long degree, long mod)
{
    if (!(0 <= val < mod) || !(1 <= degree < mod))
        perror("ValueError in is_primitive");
    long *a = unique_prime_factors(degree);
    long sz = sizeof(a)/sizeof(a[0]);


    if((power(val, degree, mod)) == 1)
    {
        for(long i=0;i<sz;i++)
        {
            if(power(val,degree/a[i],mod)==1)
                return false;
                  
        }
        return true;
    }
    return false;
    
}

long find_generator(long totient, long mod)
{
    if(!(1<=totient < mod))
        perror("Error3");
    for(long i=1;i<mod;i++)
    {
        if(is_primitive(i,totient,mod))
            return i;
    }
    perror("No generator exists");
}

long find_primitive(long degree, long totient , long mod)
{
    if(!(1 <= degree <= totient < mod))
		perror("Error1");
	if(totient % degree != 0)
		perror("Error2");
	long gen = find_generator(totient, mod);
	long root = power(gen, totient/degree,mod);

	assert(0 <= root < mod);
	return root;
}

long find_modulus(long veclen, long min)
{
    long m;
    if(veclen<1 || min<1)
        perror("ValueError in mod");
    long start = (min - 1 + veclen - 1)/veclen;
    if(start>1)
        m = start;
    else
        m=1;
    while(1)
    {
        long n = m*veclen + 1;
        assert(n>=min);
        if(is_prime(n))
            return n;
        m++;
    }
    perror("Unreachable");
}

long * transform(long* invec, long sz, long root, long mod)
{
    long c=0;
    if( (sz>=mod) || !(1<=root<mod) )
        perror("Error");
    

    long* outvec = (long*)malloc((sz+1)*sizeof(long));
    for(long i=0;i<sz;i++)
    {
        long temp=0;
        for(long j=0;j<sz;j++)
        {
            
            long val = invec[j];
            temp += ((val)* power(root,i*j,mod));
            temp %= mod;
        }
        outvec[c]=temp;
        c++;

    }
    
    return outvec;
}



long* transform_radix_2(long * invec, long sz, long root, long mod)
{
  

  long i, j, k,x;

  assert(sz > 0 && (sz & (sz - 1)) == 0); 

  j = sz>>1;
  for (i=1; i<sz-1; i++) {
    if (i < j) {
      x = invec[i]; invec[i] = invec[j]; invec[j] = x;
    }
    k = sz;
    do {
      k >>= 1;
      j ^= k;
    } while ((j & k) == 0);
  }

    long* powtable = (long*)malloc((sz+1)*sizeof(long));
    long temp=1,c=0;
    for(long i=0; i<(sz/2); i++)
    {
        powtable[c] = temp;
        temp = ((temp%mod)*(root%mod))% mod;
        c++;
    }
    long size = 2;
    while(size<=sz)
    {
        long halfsize = size/2;
        long tablestep = sz/size;
        for(long i=0; i<sz;i+=size)
        {
            long k=0;
            for(long j=i; j<(i+halfsize); j++)
            {
                long l = j+halfsize;
                long left = invec[j]%mod;
                long right = (invec[l]*powtable[k])%mod;
                invec[j] = (left + right) % mod;
				invec[l] = (left - right + mod) % mod;
				k += tablestep;
            }
        }
        size*=2;
    }
    return invec;
}

long* inverse_transform(long * invec, long sz, long root, long mod)
{
   
    long* outvec = (long*)malloc(sizeof(long)*(sz+1));
    long* out = (long*)malloc(sizeof(long)*(sz+1));
   
    outvec = transform_radix_2(invec, sz, reciprocal(root, mod), mod);
    long scaler = reciprocal(sz,mod);

    for(long i=0;i<sz;i++)
        out[i] = ((outvec[i]%mod)*(scaler%mod))% mod;
    return out;
}


long* poly_mult(long* a1, long* a2, long sz, long minmod)
{
    long* out = (long*)malloc((sz+1)*sizeof(long));
    long mod = find_modulus(sz,minmod);
    long root = find_primitive(sz,mod-1,mod);

    long* A1 = transform_radix_2(a1,sz,root,mod);
    long* A2 = transform_radix_2(a2,sz,root,mod);

    for(long i=0;i<sz;i++)
    {
        out[i] = ((A1[i])*(A2[i]))%mod;
    }
    
    return(inverse_transform(out,sz,root,mod));
   
}



int main()
{
    long n=1024,maxval=0,m1=0,m2=0,minmod,q=12289;
    long a1[1300],a2[1300];


    for(long i=0;i<n;i++)
    {
        a1[i] = rand()%(q-1);
        a2[i] = rand()%(q-1);
    }

    double time_spent = 0.0;
    clock_t begin = clock();

    for(long i=0;i<n;i++)
    {
        m1 = maxe(a1[i],m1);
        m2 = maxe(a2[i],m2);
        maxval = maxe(m1,m2); 
    }
    
    minmod = pow(maxval, 2) * n + 1;

    printf("Input vec1\n");
    for(long i=0;i<n;i++)
        printf("%ld ",a1[i]);


    printf("\nInput vec2\n");
    for(long i=0;i<n;i++)
        printf("%ld ",a2[i]);
    
    long* a = poly_mult(a1,a2,n,maxval+1);

    
    printf("\nOutput vec\n");
    for(long i=0;i<n;i++)
        printf("%ld ",a[i]);
     

    clock_t end = clock();
 

    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Elapsed time is %f seconds", time_spent);

    

}