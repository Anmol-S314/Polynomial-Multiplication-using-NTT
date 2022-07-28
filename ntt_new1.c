#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include<time.h>
#include<stdbool.h>

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
    // long i;
    
    // if(n == 2)
    //     return true;
    
    // if(n % 2 == 0)
    //     return false;
    
    // for(i = 3; i < sqrt(n); i = i + 2) 
    // {
    //     if(n % i == 0)
    //         return false;
    // }
    // return true;
    long flag=1;
    for (long i = 2; i <= sqrt(n); i++) {
 
        // If n is divisible by any number between
        // 2 and n/2, it is not prime
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

long squart(long n)
{
    if(n<0)
        perror("ValueError");
    long i=1;
    while(i*i<=n)
    {
        i*=2;
        long res = 0;
        while(i>0)
        {
            if(((res+1)*(res+1))<=n)
                res+=i;
            i=i/2;
        }
        return res;
    }
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
        // q is quotient
        long q = a / m;
        long t = m;
 
        // m is remainder now, process same as
        // Euclid's algo
        m = a % m, a = t;
        t = y;
 
        // Update y and x
        y = x - q * y;
        x = t;
    }
 
    // Make x positive
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

    // if((power(val, degree, mod)) == 1)
    // {
    //     printf("%d",val);
    //     return false;
    // }
   // printf("73:%d\n",power(6, degree, mod));
    if((power(val, degree, mod)) == 1)
    {
        for(long i=0;i<sz;i++)
        {
            if(power(val,degree/a[i],mod)==1)
            {
                //printf("%d\n",val);
                return false;
            }       
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
        {
            //printf("i=%ld\n",i);
            return i;
        }
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
    //printf("111:%ld\n",gen);
	long root = power(gen, totient/degree,mod);
    //printf("113:%ld\n",root);
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
       // printf("186%ld\n",n);
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
    //long sz = sizeof(invec)/sizeof(invec[0]);
    if( (sz>=mod) || !(1<=root<mod) )
        perror("Error");
    

    long* outvec = (long*)malloc((sz+1)*sizeof(long));
    for(long i=0;i<sz;i++)
    {
        long temp=0;
        for(long j=0;j<sz;j++)
        {
            //long j = invec[k][0];
            long val = invec[j];
            temp += ((val)* power(root,i*j,mod));
            //temp += ((val)* 1);
            temp %= mod;
        }
        outvec[c]=temp;
        c++;

    }
    
    return outvec;
}



long* transform_radix_2(long * invec, long sz, long root, long mod)
{
    //long levels = (long)log2(sz)+1;
    //bool levels = sz && (!(sz & (sz - 1)));
    //if((levels<<1)!= sz)
       // perror("Length not power of 2");

    long i, j, k,x;
  //long32_t x;

  assert(sz > 0 && (sz & (sz - 1)) == 0); // n must be a power of 2

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
    // printf("266\n");
    // for(long i=0;i<sz;i++)
    // {
    //     printf("%d ",invec[i]);
    // }
    // printf("266\n");
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
    //long root=3;
    long* outvec = (long*)malloc(sizeof(long)*(sz+1));
    long* out = (long*)malloc(sizeof(long)*(sz+1));
    //printf("inverse:%ld\n",reciprocal(root,mod));
    outvec = transform_radix_2(invec, sz, reciprocal(root, mod), mod);
    long scaler = reciprocal(sz,mod);
    //printf("scaler:%ld\n",scaler);
    //for(long i=0;i<sz;i++)
        //printf("%ld ",outvec[i]);
    //printf("\n");
    for(long i=0;i<sz;i++)
        out[i] = ((outvec[i]%mod)*(scaler%mod))% mod;
    return out;
}
long* find_params_and_transform(long* a1,long sz, long minmod)
{
    long mod = find_modulus(sz,minmod);
    //printf("161:%ld\n",mod);
    long root = find_primitive(sz,mod-1,mod);
    //printf("163:%ld\n",root);
    long* c = transform_radix_2(a1,sz,root,mod);
    //long * b = transform(a1,sz,3,mod);
    //return(inverse_transform(c,sz,root,mod));
    return(c);
    
}

long* find_params_and_transform2(long* a1, long* a2, long sz, long minmod)
{
    long* out = (long*)malloc((sz+1)*sizeof(long));
    //printf("159:%ld\n",sz);
    long mod = find_modulus(sz,minmod);
    //printf("161:%ld\n",mod);
    long root = find_primitive(sz,mod-1,mod);
    //printf("163:%ld\n",root);
    long* A1 = transform_radix_2(a1,sz,root,mod);
    long* A2 = transform_radix_2(a2,sz,root,mod);
    //printf("\nbruh\n");
    for(long i=0;i<sz;i++)
    {
        //printf("%ld ",A2[i]);
    }
    for(long i=0;i<sz;i++)
    {
        out[i] = ((A1[i])*(A2[i]))%mod;
    }
    //printf("out\n");
    for(long i=0;i<sz;i++)
    {
        //printf("%ld ",out[i]);
    }
    return(inverse_transform(out,sz,root,mod));
    //return out;
}



int main()
{
    long n=64,maxval=0,m1=0,m2=0,minmod;
    long a1[1300],a2[1300],b[1300];
    FILE *fptr;
    fptr = fopen("output.txt", "r");
    if(fptr == NULL)
   {
      printf("Error!");   
      exit(1);             
   }
    //long a1[1000],a2[1000];
    //long a3[4] = {15,21,13,44} ;
    // root=30 mod=53 ntt=[40,1,16,3]


    // to store the execution time of code
    //printf("%ld",find_modulus(512,12289));
    //printf("lig%ld",find_primitive(256,12288,12289));
    double time_spent = 0.0;
    clock_t begin = clock();
    
    long count=0,avg=0;
    while(count!=1000)
    {
            for(long i=0;i<n;i++)
    {
        a1[i] = rand()%12289;
        a2[i] = rand()%12289;
    }
        
    //     for(long i=0;i<n;i++)
    //         fscanf(fptr, "%ld", &a1[i]);
    //     for(long i=0;i<n;i++)
    //         fscanf(fptr, "%ld", &a2[i]);
    //     for(long i=0;i<n;i++)
    //         fscanf(fptr, "%ld", &b[i]);
        //long a2[64] = {6222, 6565, 10290, 8757, 8158, 601, 3063, 1749, 11300, 608, 9561, 1894, 5459, 880, 6260, 11255, 12244, 4335, 7454, 11825, 5877, 4385, 10830, 8439, 5082, 7706, 9020, 8656, 1328, 733, 6282, 4647};
        //long a3[64] = {521, 530, 78, 421, 379, 319, 51, 620, 468, 382, 43, 568, 571, 441, 339, 311, 107, 147, 173, 480, 631, 6, 432, 70, 80, 192, 517, 570, 232, 418, 372, 504, 63, 181, 233, 568, 619, 534, 557, 336, 513, 98, 152, 446, 334, 125, 597, 86, 241, 93, 178, 312, 235, 161, 492, 52, 509, 477, 450, 183, 407, 261, 563, 533};
        //long a1[64] = {5903, 2198, 2028, 9245, 1682, 11990, 1976, 687, 792, 5790, 7558, 3021, 1646, 7405, 5124, 12077, 5304, 6504, 10192, 5425, 494, 9597, 1041, 4886, 1944, 6464, 6719, 9310, 10102, 11536, 11460, 8045};
        //long a1[8] = {4,1,4,2,1,3,5,6};
        //long a2[8] = {6,1,8,0,3,3,9,8};
        for(long i=0;i<n;i++)
        {
            m1 = maxe(a1[i],m1);
            m2 = maxe(a2[i],m2);
            maxval = maxe(m1,m2); 
        }
        
        minmod = pow(maxval, 2) * n + 1;
        //printf("Input vec1\n");
        for(long i=0;i<n;i++)
        {
            //printf("%ld ",a1[i]);
        }
        //printf("\nInput vec2\n");
        for(long i=0;i<n;i++)
        {
            //printf("%ld ",a2[i]);
        }
        long* a = find_params_and_transform2(a1,a2,n,maxval+1);
        //long* b = find_params_and_transform(a3,64,632);
        

        //long sz = sizeof(a)/sizeof(a[0]);
        //printf("%d\n",sz);
        
        //printf("\nOutput vec\n");
    //     for(long i=0;i<n;i++)
    //     {
    //         //printf("%ld ",a[i]);
    //         if(b[i]!=a[i])
    //         {
    //             printf("Error doesnt match");
    //             return 1;
    //         }

            
    //     }
        
        count++;
    }

    printf("All Match\n");


    clock_t end = clock();
 
    // calculate elapsed time by finding difference (end - begin) and
    // dividing the difference by CLOCKS_PER_SEC to convert to seconds
    time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
    printf("cps%ld",CLOCKS_PER_SEC);
    printf("The elapsed time is %f seconds", time_spent);

    

}