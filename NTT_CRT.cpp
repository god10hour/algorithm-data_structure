#include <bits/stdc++.h>
using namespace std;

template<int mod, int root>
struct NTT {
    static int pow(int n, int k) {
        int r=1;
        while (k) {
            if (k&1)
                r=1LL*r*n%mod;
            n=1LL*n*n%mod;
            k/=2;
        }
        return r;
    }

    static void ntt(vector<int> &v, bool inv) {
        const int s=v.size();
        for (int i=1, j=0; i<s; i++) {
            int bit=s>>1;
            while (!((j^=bit)&bit)) bit>>=1;
            if (i<j) swap(v[i], v[j]);
        }
        for (int k=1; k<s; k*=2) {
            int w=pow(root, (mod-1)/(k*2));
            if (inv) w=pow(w, mod-2);
            for (int i=0; i<s; i+=k*2) {
                int uni=1;
                for (int j=0; j<k; j++) {
                    int a=1LL*v[i+j];
                    int b=1LL*v[i+j+k]*uni%mod;
                    v[i+j]=(1LL*a+b)%mod;
                    v[i+j+k]=(1LL*a-b+mod)%mod;
                    uni=1LL*uni*w%mod;
                }
            }
        }
        if (inv) {
            int I=pow(s, mod-2);
            for (int i=0; i<s; i++)
                v[i]=1LL*v[i]*I%mod;
        }
    }

    static vector<int> conv(vector<int> v, vector<int> u) {
        if (v.empty() || u.empty()) return {};
        int s=1, n=v.size()+u.size()-1;
        while (s<n) s*=2;
        v.resize(s); ntt(v, 0);
        u.resize(s); ntt(u, 0);
        for (int i=0; i<s; i++)
            v[i]=1LL*v[i]*u[i]%mod;
        ntt(v, 1), v.resize(n);
        return v;
    }
};

struct CRT {
    using u128=__uint128_t;
    static constexpr int MODS[3]={
        998244353, 469762049, 1224736769
    };
    static constexpr int ROOT=3;
    NTT<MODS[0], ROOT> P0;
    NTT<MODS[1], ROOT> P1;
    NTT<MODS[2], ROOT> P2;

    static int pow(int n, int k, int mod) {
        int r=1;
        while (k) {
            if (k&1)
                r=1LL*r*n%mod;
            n=1LL*n*n%mod;
            k/=2;
        }
        return r;
    }

    int64_t combine(const int a[3]) {
        u128 r=0, mod=1;
        for (int i=0; i<3; i++) {
            int64_t Mi=MODS[i], ri=r%Mi;
            if (ri<0) ri+=Mi;
            int diff=(1LL*a[i]+Mi-ri)%Mi;
            int inv=pow(mod%Mi, Mi-2, Mi);
            r+=mod*(diff*inv%Mi);
            mod*=Mi;
        }
        return r;
    }

    vector<int> conv(const vector<int> &V, const vector<int> &U, int p) {
        if (V.empty() || U.empty()) return {};
        int n=V.size()+U.size()-1;
        array<vector<int>, 3> v, u;
        for (int i=0; i<3; i++) v[i].resize(V.size());
        for (int i=0; i<3; i++) u[i].resize(U.size());

        for (int i=0; i<V.size(); i++)
            for (int j=0; j<3; j++) {
                v[j][i]=V[i]%MODS[j];
                if (v[j][i]<0) v[j][i]+=MODS[j];
            }
        for (int i=0; i<U.size(); i++)
            for (int j=0; j<3; j++) {
                u[j][i]=U[i]%MODS[j];
                if (u[j][i]<0) u[j][i]+=MODS[j];
            }

        auto w0=P0.conv(v[0], u[0]);
        auto w1=P1.conv(v[1], u[1]);
        auto w2=P2.conv(v[2], u[2]);

        vector<int> r(n);
        for (int i=0; i<n; i++) {
            int a[3]={
                i<w0.size()? w0[i]:0,
                i<w1.size()? w1[i]:0,
                i<w2.size()? w2[i]:0
            };
            r[i]=combine(a)%p;
        }
        return r;
    }
};

#define int int64_t
template<int mod=0, int root=0>
struct Poly {
    inline static int lim=-1;
    vector<int> v;
    Poly()=default;
    explicit Poly(int n): v(n, 0) {}
    explicit Poly(const vector<int> &v): v(v) {}
    explicit Poly(vector<int> &&v): v(move(v)) {}

    [[nodiscard]] size_t size() const { return v.size(); }
    void resize(int n) { v.resize(n); }
    int &operator[](int i) { return v[i]; }
    void selfcut(int m) { if (m<(int)v.size()) v.resize(m); }
    vector<int> cut(int m) {
        if (m>=(int)v.size()) return v;
        return vector<int>(v.begin(), v.begin()+m);
    }
    static void set(int m) { lim=m; }
    static void clear() { lim=-1; }
    void trim() {
        while (!v.empty() && v.back()==0)
            v.pop_back();
    }

    static Poly conv(const Poly &V, const Poly &U) {
        Poly r;
        if constexpr (mod)
            r.v=NTT<mod, root>::conv(V.v, U.v);
        if (lim>0 && r.v.size()>lim)
            r.resize(lim);
        return r;
    }

    Poly operator+(const Poly &o) const {
        int s=max(size(), o.size());
        if (lim>0) s=min(s, lim);
        Poly r(s);
        for (int i=0; i<s; i++) {
            int x=i<size()? v[i]:0;
            int y=i<o.size()? o.v[i]:0;
            if constexpr (mod) r[i]=(x+y)%mod;
            else r[i]=x+y;
        }
        return r;
    }
    Poly operator-(const Poly &o) const {
        int s=max(size(), o.size());
        if (lim>0) s=min(s, lim);
        Poly r(s);
        for (int i=0; i<s; i++) {
            int x=i<size()? v[i]:0;
            int y=i<o.size()? o.v[i]:0;
            if constexpr (mod) r[i]=(x-y+mod)%mod;
            else r[i]=x-y;
        }
        return r;
    }
    Poly operator*(const Poly &o) const {
        if (lim>0) {
            Poly w=conv(*this, o);
            if (lim>0 && w.size()>lim) w.resize(lim);
            return w;
        }
        return conv(*this, o);
    }

    Poly pow(Poly base, int k) {
        Poly r(1);
        r[0]=1;
        while (k) {
            if (k&1)
                r=r*base;
            base=base*base;
            k/=2;
        }
        if (lim>0 && base.size()>lim)
            base.resize(lim);
        return r;
    }
};
#undef int
