#include <bits/stdc++.h>
using namespace std;

/*  example
 *    Nimber64 a(2), b(3);
 *    auto c = a * b;       // nimber multiplication
 *    auto s = a + b;       // XOR
 *    auto inv = a.inv();   // inverse
 */

class Nimber {
public:
    using u16=uint16_t;
    using u32=uint32_t;
    using u64=uint64_t;
    using u128=__uint128_t;

private:
    static constexpr u16 G16=10279U;

    static inline u16 expBuffer_[4*(1u<<16)+16];
    static inline u16 *exp_=nullptr;
    static inline u16 *exp3=nullptr;
    static inline u16 *exp6=nullptr;
    static inline int log_[1U<<16];

    static inline bool inited=false;
    u64 v=0;

    template<int L>
    static u64 star_slow(u64 a, u64 b) {
        static_assert(L>0);
        if constexpr (L==1) return a&b;
        else {
            constexpr int l=L>>1;
            const u64 mask=(1ULL<<l)-1;
            const u64 a0=a&mask, a1=a>>l;
            const u64 b0=b&mask, b1=b>>l;

            const u64 ab=star_slow<l>(a0, b0);
            const u64 lo=ab^star_slow<l>(1ULL<<(l-1), star_slow<l>(a1, b1));
            const u64 hi=ab^star_slow<l>(a0^a1, b0^b1);
            return lo | hi<<l;
        }
    }

    static u16 star16(u16 a, u16 b) {
        if (!a || !b) return 0;
        return exp_[log_[a]+log_[b]];
    }

    static u16 star16offset(u16 a, u16 b, int offset) {
        if (!a || !b) return 0;
        return exp_[log_[a]+log_[b]+offset];
    }

    static u32 star32(u32 a, u32 b) {
        const u16 a0=u16(a), a1=u16(a>>16);
        const u16 b0=u16(b), b1=u16(b>>16);
        const u16 a01=a0^a1, b01=b0^b1;

        const u16 ab=star16(a0, b0);
        const u16 lo=ab^star16offset(a1, b1, 3);
        const u16 hi=ab^star16(a01, b01);
        return u32(lo) | u32(hi)<<16;
    }

    static u32 omega(u32 a) {
        const u16 a0=u16(a), a1=u16(a>>16);
        const u16 a01=u16(a0^a1);

        const u16 lo=a1? exp6[log_[a1]]: 0;
        const u16 hi=a01? exp3[log_[a01]]: 0;
        return u32(lo) | u32(hi)<<16;
    }

    static u64 star64(u64 a, u64 b) {
        const u32 a0=u32(a), a1=u32(a>>32);
        const u32 b0=u32(b), b1=u32(b>>32);
        const u32 a01=a0^a1, b01=b0^b1;

        const u32 a0b0=star32(a0, b0);
        const u32 lo=a0b0^omega(star32(a1, b1));
        const u32 hi=a0b0^star32(a01, b01);
        return u64(lo) | u64(hi)<<32;
    }

    static void init() {
        if (inited) return;

        u16 *center=expBuffer_+(2*(1u<<16)+8);
        exp_=center;
        exp3=exp_+3;
        exp6=exp_+6;

        exp_[0]=1;
        constexpr int PERIOD=(1<<16)-1;
        for (int i=1; i<PERIOD; i++)
            exp_[i]=star_slow<16>(exp_[i-1], G16);
        for (int i=PERIOD; i<2*(1<<16); i++)
            exp_[i]=exp_[i-PERIOD];

        for (int i=0; i<PERIOD; i++)
            log_[exp_[i]]=i;
        log_[0]=-INT_MAX/2;

        inited=true;
    }

public:
    Nimber(): v(0) { init(); }
    explicit Nimber(u64 x): v(x) { init(); }

    u64 value() const { return v; }
    bool isZero() const { return v==0; }

    friend Nimber operator^(Nimber a, Nimber b) { return Nimber(a.v^b.v); }
    friend Nimber operator+(Nimber a, Nimber b) { return a^b; }
    friend Nimber operator-(Nimber a, Nimber b) { return a^b; }

    Nimber &operator^=(Nimber o) { v^=o.v; return *this; }
    Nimber &operator+=(Nimber o) { v^=o.v; return *this; }
    Nimber &operator-=(Nimber o) { v^=o.v; return *this; }

    friend Nimber operator*(Nimber a, Nimber b) {
        return Nimber(star64(a.v, b.v));
    }
    Nimber &operator*=(Nimber o) { v=star64(v, v*o.v); return *this; }

    bool operator==(Nimber o) const { return v==o.v; }
    bool operator!=(Nimber o) const { return v!=o.v; }

    Nimber pow(u128 k) const {
        Nimber base=*this, r=Nimber(1);
        while (k) {
            if (k&1)
                r*=base;
            base*=base;
            k>>=1;
        }
        return r;
    }

    Nimber inv() const {
        if (isZero())
            throw runtime_error("Nimber64: inverse of zero");
        u128 exp=((u128)1<<64)-2;
        return pow(exp);
    }

    friend Nimber operator/(Nimber a, Nimber b) {
        if (b.isZero())
            throw runtime_error("Nimber64: divide by zero");
        return a*b.inv();
    }
    Nimber &operator/=(Nimber b) {
        *this=*this/b; return *this;
    }
};
