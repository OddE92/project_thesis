#ifndef RAN_GEN
#define RAN_GEN

/***** DISCLAIMER *****/
/*      This is the random number generator given in:
    Numerical Recipes, third edition, Cambridge University Press 2007
        Press, William H.
        Teukolsky, Saul A.
        Vetterling, William T.
    p. 342-343

    It sets up a struct, where the object is initialized with a seed. Basic use is as follows:

    Ran object(SEED);

    object.int64() returns a 64-bit integer
    object.int32() returns a 32-bit integer
    object.doub()  returns a double in the range [0, 1]

    object.doub(); object.doub() returns two different numbers consecutively.

    All credit is given to the authors.
*/

typedef unsigned long long int Ullong;
typedef unsigned int Uint;

struct Ran
{

    Ullong u, v, w;
    Ran(Ullong j) : v(4101842887655102017LL), w(1)                      //Constructor.Call with any integer seed(except value of v above).
    { 
        u = j ^ v;
        int64();
        v = u;
        int64();
        w = v;
        int64();
    }
    inline Ullong int64()                                               //Return 64-bit random integer.
    {
        u = u * 2862933555777941757LL + 7046029254386353087LL;
        v ^= v >> 17;
        v ^= v << 31;
        v ^= v >> 8;
        w = 4294957665U * (w & 0xffffffff) + (w >> 32);
        Ullong x = u ^ (u << 21);
        x ^= x >> 35;
        x ^= x << 4;
        return (x + v) ^ w;
    }
    inline double doub() { return 5.42101086242752217E-20 * int64(); }  //Return random double-precision floating value in the range 0. to 1.
    inline Uint int32() { return (Uint)int64(); }                       //Return random 32-bit integer
};

#endif