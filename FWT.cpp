#include <bits/stdc++.h>
using namespace std;

enum class mode { XOR, OR, AND };
template<mode OP, int mod>
struct FWT {
    inline static int inv2=(mod+1)/2;
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

    static void fwt(vector<int> &v, bool inv) {
        const int s=v.size();
        for (int k=1; k<s; k*=2)
            for (int i=0; i<s; i+=k*2)
                for (int j=0; j<k; j++) {
                    int a=v[i+j], b=v[i+j+k];
                    if constexpr (OP==mode::XOR) {
                        v[i+j]=(a+b)%mod;
                        v[i+j+k]=(a-b+mod)%mod;
                    }
                    if constexpr (OP==mode::OR) {
                        if (!inv)
                            v[i+j+k]=(a+b)%mod;
                        else v[i+j+k]=(b-a+mod)%mod;
                    }
                    if constexpr (OP==mode::AND) {
                        if (!inv)
                            v[i+j]=(a+b)%mod;
                        else
                            v[i+j]=(a-b+mod)%mod;
                    }
                }
        if constexpr (OP==mode::XOR)
            if (inv) {
                int levels=0;
                for (int i=s; i>1; i/=2) levels++;
                int I=pow(inv2, levels);
                for (int i=0; i<s; i++)
                    v[i]=1LL*v[i]*I%mod;
            }
    }

    static vector<int> conv(vector<int> v, vector<int> u) {
        if (v.empty() || u.empty()) return {};
        int s=1, n=max(v.size(), u.size());
        while (s<n) s*=2;
        vector<int> w(s);
        v.resize(s), fwt(v, 0);
        u.resize(s), fwt(u, 0);
        for (int i=0; i<s; i++)
            w[i]=1LL*v[i]*u[i]%mod;
        fwt(w, 1);
        return w;
    }
};
