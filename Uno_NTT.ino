#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include<time.h>
#include<stdbool.h>

/* The code performs polynomial multiplication for array sizes upto s=64 with the modulus at q=12289 . 
 * For faster computation the root values are precomputed in an array but this could
 * also be modified to compute any mod and root value given an array value.
 * To use the code simply plug in the array size which must be a power of 2 and it will compute the given multiplication.
 */
template <typename T>
Print& operator<<(Print& printer, T value)
{
  printer.print(value);
  return printer;
}


long maxe(long x, long y)
{
  if (x >= y)
    return x;
  else
    return y;
}


long power(long x, long y, long mod)
{
  long z = 1;
  for (long i = 0; i < y; i++)
    z = ((z % mod) * (x % mod)) % mod;
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


long* unique_prime_factors(long n)
{
  long* prime = new long(20);
  long sz = 0, i = 2, end = sqrt(n);

  if (n < 1)
    return NULL;
  while (i <= end)
  {
    if (n % i == 0)
    {
      n = n / i;
      prime[sz] = i;
      sz++;
      while (n % i == 0)
        n = n / i;
      end = sqrt(n);
    }
    i += 1;
  }
  if (n > 1)
  {
    prime[sz] = n;
    sz++;
  }

  prime[sz] = 0;
  sz++;
  
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

  // Make x positive
  if (x < 0)
    x += m0;

  return x;
}


bool is_primitive(long val, long degree, long mod)
{
  long i = 0;
  if (!(0 <= val < mod) || !(1 <= degree < mod))
  {
    Serial << "Er1";
    exit(0);
  }
  long* prime = unique_prime_factors(degree);


  if ((power(val, degree, mod)) == 1)
  {
    while (prime[i] != 0)
    {
      if (power(val, degree / prime[i], mod) == 1)
        return false;
      i++;
    }
    return true;
  }
  return false;

}


long find_generator(long totient, long mod)
{
  if (!(1 <= totient < mod))
  {
    Serial << "Er2";
    exit(0);
  }
  for (long i = 2; i < mod; i++)
  {
    if (is_primitive(i, totient, mod))
      return i;
  }
  Serial << "Er3";
  exit(1);
}

long find_primitive(long degree, long totient , long mod)
{
  if (!(1 <= degree <= totient < mod))
  {
    Serial << "Er4";
    exit(0);
  }
  if (totient % degree != 0)
  {
    Serial << "Er5";
    exit(1);
  }
  long gen = find_generator(totient, mod);
  long root = power(gen, totient / degree, mod);
  
  assert(0 <= root < mod);
  return root;
}

long find_modulus(long veclen, long min)
{
  long m;
  if (veclen < 1 || min < 1)
  {
    Serial << "Er6";
    exit(0);
  }
  long start = (min - 1 + veclen - 1) / veclen;
  if (start > 1)
    m = start;
  else
    m = 1;
  for (long i = m;; i++)
  {
    long n = i * veclen + 1;
    assert(n >= min);
    if (is_prime(n))
      return n;
  }
  Serial << "Er7";
  exit(1);
}


long * transform(long* invec, long sz, long root, long mod)
{
    long c=0;
    if( (sz>=mod) || !(1<=root<mod) )
        perror("Error");


    long* outvec = new long[100];
    for(long i=0;i<sz;i++)
    {
        long temp=0;
        for(long j=0;j<sz;j++)
        {
            long val = invec[j];
            temp += ((val) * power(root,i*j,mod));
            temp %= mod;
        }
        outvec[c]=temp;
        c++;

    }

    return outvec;
}




long* transform_radix_2(long* invec, long sz, long root, long mod)
{
  long i, j, k, x;
  assert(sz > 0 && (sz & (sz - 1)) == 0); // n must be a power of 2

  j = sz >> 1;
  for (i = 1; i < sz - 1; i++) {
    if (i < j) {
      x = invec[i]; invec[i] = invec[j]; invec[j] = x;
    }
    k = sz;
    do {
      k >>= 1;
      j ^= k;
    } while ((j & k) == 0);
  }

  long* powtable = new long[100];
  long temp = 1, c = 0;
  for (long i = 0; i < (sz / 2); i++)
  {
    powtable[c] = temp;
    temp = ((temp) * (root)) % mod;
    c++;
  }
 
  long size = 2;
  while (size <= sz)
  {
    long halfsize = size / 2;
    long tablestep = sz / size;
    for (long i = 0; i < sz; i += size)
    {
      long k = 0;
      for (long j = i; j < (i + halfsize); j++)
      {
        long l = j + halfsize;
        long left = invec[j] % mod;
        long right = (invec[l] * powtable[k]) % mod;
        invec[j] = (left + right) % mod;
        invec[l] = (left - right + mod) % mod;
        k += tablestep;
      }
    }
    size *= 2;
  }
  delete [] powtable; 
  return invec;
}


long* find_params_and_transform2(long* a1, long* a2, long sz, long minmod)
{
  long* out = new long[140];
  //long* mem = new mem[10];
  //long* A1 = new long[100];
  //long* A2 = new long[100];

  //long mod = find_modulus(sz, minmod);
  
  long mod=12289;
  
  //long root=5860;
  //long root=7311;
  //long root=12149;
  //long root=8340;
  
  //long root = find_primitive(sz, mod - 1, mod);
  long mem[9] = {1,12288,1479,8246,4134,5860,7311,12149,8340};
  //long mem[9] = {1,768,707,173,712,251,64,243,8340};
  long temp = log(sz)/log(2);
  long root = mem[temp];

  Serial<<"mod,root"<<mod<<root<<"\n";

  long int t1 = millis();
  long* A1 = transform_radix_2(a1, sz, root, mod);
  long* A2 = transform_radix_2(a2, sz, root, mod);

  for (long i = 0; i < sz; i++)
  {
    out[i] = ((A1[i]) * (A2[i]))%mod;
    
  }
  
  //delete [] A1;
  //delete [] A2;
  long tem = reciprocal(root, mod);
  //long tem = 8747;
  //long tem = 9650;
  out = transform_radix_2(out,sz,tem, mod);
  long scaler = reciprocal(sz, mod);
  //long scaler=11905;
  //long scaler=12097;
  
  Serial<<tem<<" "<<scaler;
  
  for (long i = 0; i < sz; i++)
    out[i] = ((out[i] % mod) * (scaler % mod)) % mod;
   long int t2 = millis();

  Serial.println("NTT Time:");
  Serial.println((t2 - t1));
    
  return out;
  
}


void setup()
{

  Serial.begin(9600);
  long int t1 = millis();
  long n=64,maxval = 0, m1 = 0, m2 = 0, minmod,q=12289,count=0;
  long* a1 = new long[64];
  long* a2 = new long[64];
  randomSeed(analogRead(0));


    for(long i=0;i<n;i++)
    {
        a1[i] = random()%(q-1);
        a2[i] = random()%(q-1);
        //Serial<<"a,b"<<a<<b<<"\n";
        m1 = maxe(a1[i], m1);
        m2 = maxe(a2[i], m2);
        maxval = maxe(m1, m2);
    }
  
  
  
    Serial.print("Input vec1\n");
    for (long i = 0; i < n; i++)
    {
      Serial << a1[i] << ", ";
    }
    
    Serial.print("\nInput vec2\n");
    for (long i = 0; i < n; i++)
    {
      Serial << a2[i] << ", ";
    }
    
    long* a = find_params_and_transform2(a1, a2, n, maxval+1);
    
  
    Serial.print("\nOutput vec\n");
    for (long i = 0; i < n; i++)
    {
      Serial << a[i] << ", "; 
    }
    Serial.print("\n");
    
  
  long int t2 = millis();
  delete [] a1;
  delete [] a2;
  Serial.println("Time:");
  Serial.println((t2 - t1));


}


void loop(){}
