#include <bits/stdc++.h>
using namespace std;
using ldb=long double;

template<typename Type>
struct Vec {
    Type x, y, z;
    Vec()=default;
    explicit Vec(Type x, Type y, Type z): x(x), y(y), z(z) {}

    Vec operator+(const Vec &o) const { return Vec(x+o.x, y+o.y, z+o.z); }
    Vec operator-(const Vec &o) const { return Vec(x-o.x, y-o.y, z-o.z); }
    Vec operator*(Type k) const { return Vec(x*k, y*k, z*k); }
    Vec& operator+=(const Vec& o) { x+=o.x; y+=o.y; z+=o.z; return *this; }

    friend Type operator*(const Vec &a, const Vec &b) {
        return a.x*b.x + a.y*b.y + a.z*b.z;
    }
    friend  Vec operator/(const Vec &a, const Vec &b) {
        return Vec(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
    }

    static Type norm2(const Vec &a) { return a*a; }
    static Type norm(const Vec &a) { return sqrtl(norm2(a)); }

    static Vec normalize(const Vec &a) {
        Type n=norm(a);
        return n==0? a: a*(1.0L/n);
    }
};

template<typename Type>
struct Convex3D {
    using V=Vec<Type>;
    struct Face { int a, b, c; };
    vector<Vec<Type>> P;
    vector<Face> faces;
    const Type EPS;
    Convex3D(Type eps=1e-12L): EPS(eps) {}
    explicit Convex3D(vector<Vec<Type>> pts, Type eps=1e-12L): EPS(eps) { get(move(pts)); }

    static Type s_dist(const V &p, const V &A, const V &B, const V &C) {
        Vec n=(B-A)/(C-A);
        return n*(p-A);
    }

    bool init(vector<int> &idx) {
        int n=P.size();
        int i0=-1, i1=-1, i2=-1, i3=-1;
        for (int i=0; i<n; i++) {
            if (i0==-1) i0=i;
            else if (V::norm(P[i]-P[i0])>EPS) {
                i1=i;
                break;
            }
        }
        if (i1==-1) return false; // all same
        for (int i=i1+1; i<n; i++) {
            V u=P[i1]-P[i0], v=P[i]-P[i0];
            if (V::norm(u/v)>EPS) {
                i2=i;
                break;
            }
        }
        if (i2==-1) return false; // collinear
        for (int i=i2+1; i<n; i++)
            if (fabsl(s_dist(P[i], P[i0], P[i1], P[i2]))>EPS) {
                i3=i;
                break;
            }
        if (i3==-1) return false; // coplanar
        idx={i0, i1, i2, i3};
        return true;
    }

    void normalize_face(Face &f, int q) {
        const V &A=P[f.a], &B=P[f.b], &C=P[f.c];
        if (s_dist(P[q], A, B, C)>0)
            swap(f.b, f.c);
    }

    void get(vector<V> pts) {
        P=move(pts);
        faces.clear();
        int n=P.size();
        if (n<4) return;

        vector<int> idx;
        if (!init(idx)) return;

        int i0=idx[0], i1=idx[1], i2=idx[2], i3=idx[3];
        vector<Face> init={
            {i0, i1, i2},
            {i0, i3, i1},
            {i0, i2, i3},
            {i1, i3, i2}
        };
        for (auto &f: init)
            normalize_face(f, i3);
        faces=init;

        for (int p=0; p<n; p++) {
            bool any=false;
            for (auto &f: faces)
                if (s_dist(P[p], P[f.a], P[f.b], P[f.c])>EPS) {
                    any=true;
                    break;
                }
            if (!any) continue;

            vector<int> vis; vis.reserve(faces.size());
            for (int i=0; i<faces.size(); i++) {
                auto &f=faces[i];
                if (s_dist(P[p], P[f.a], P[f.b], P[f.c])>EPS)
                    vis.push_back(i);
            }

            map<pair<int, int>, int> cnt;
            vector<bool> use(faces.size(), 1);
            for (int i: vis) use[i]=0;
            auto push=[&](int u, int v) {
                cnt[{u, v}]++;
            };

            for (int i: vis) {
                auto &f=faces[i];
                push(f.a, f.b);
                push(f.b, f.c);
                push(f.c, f.a);
            }

            vector<pair<int, int>> horizon; horizon.reserve(cnt.size());
            for (auto &kv: cnt) {
                auto e=kv.first;
                auto re=make_pair(e.second, e.first);
                int c=kv.second;
                int rc=cnt.count(re)? cnt[re]: 0;
                if (c==1 && rc==0) horizon.push_back(e);
            }

            vector<Face> nf; nf.reserve(faces.size());
            for (int i=0; i<faces.size(); i++) if (use[i])
                nf.push_back(faces[i]);
            for (auto &e: horizon) {
                Face add{e.first, e.second, p};
                int q=e.first;
                const Vec<Type> &A=P[add.a], &B=P[add.b], &C=P[add.c];
                if (s_dist(P[q], A, B, C)>0)
                    swap(add.b, add.c);
                nf.push_back(add);
            }
            faces.swap(nf);
        }

        auto key=[&](Face f) {
            array<int, 3> s{f.a, f.b, f.c};
            sort(s.begin(), s.end());
            return s;
        };
        set<array<int, 3>> seen;
        vector<Face> uniq; uniq.reserve(faces.size());
        for (auto &f: faces)
            if (seen.insert(key(f)).second)
                uniq.push_back(f);
        faces.swap(uniq);

        V G(0, 0, 0);
        for (auto &v: P) G+=v;
        if (!P.empty()) G=G*(Type(1)/P.size());
        for (auto &f: faces) {
            const Vec<Type> &A=P[f.a], &B=P[f.b], &C=P[f.c];
            if (s_dist(G, A, B, C)>0)
                swap(f.b, f.c);
        }
    }

    ldb surface() const {
        ldb S=0;
        for (auto &f: faces) {
            Vec<Type> a=P[f.a], b=P[f.b], c=P[f.c];
            S+=0.5L*V::norm((b-a)/(c-a));
        }
        return S;
    }

    ldb volume() const {
        ldb s=0;
        for (auto &f: faces) {
            Vec<Type> a=P[f.a], b=P[f.b], c=P[f.c];
            s+=a*(b/c); // det(a, b, c)
        }
        return fabsl(s)/6.0L;
    }

    vector<int> indices() const {
        vector<bool> use(P.size(), 0);
        vector<int> idx; idx.reserve(P.size());
        for (auto &f: faces)
            use[f.a]=use[f.b]=use[f.c]=1;
        for (int i=0; i<use.size(); i++) if (use[i])
            idx.push_back(i);
        return idx;
    }
};

int main() {
    cin.tie(0)->sync_with_stdio(0);
    int n;
    cin >> n;
    vector<Vec<ldb>> P(n);
    for (int i=0; i<n; i++) {
        ldb x, y, z;
        cin >> x >> y >> z;
        P[i]=Vec(x, y, z);
    }
    Convex3D H(P);
    auto idx=H.indices();

    cout << idx.size() << '\n';
    for (int x: idx) cout << x << ' ';
    cout << '\n';

    cout << H.volume() << '\n';
    cout << H.surface() << '\n';
}
