#include <bits/stdc++.h>
using namespace std;

template <typename Type>
struct Poi {
    Type x, y;
    Poi()=default;
    explicit Poi(Type x, Type y): x(x), y(y) {}

    Poi operator+(const Poi &p) const { return Poi(x+p.x, y+p.y); }
    Poi operator-(const Poi &p) const { return Poi(x-p.x, y-p.y); }
    Poi operator*(Type k)       const { return Poi(x*k, y*k); }
    friend Poi operator*(Type k, const Poi &p) { return Poi(p.x*k, p.y*k); }

    Type operator*(const Poi &p) const { return x*p.x + y*p.y; }
    Type operator/(const Poi &p) const { return x*p.y - y*p.x; }

    bool operator==(const Poi &p) const { return x==p.x && y==p.y; }
    bool operator< (const Poi &p) const { return x<p.x || (x==p.x && y<p.y); }

    static Type dis(const Poi &a, const Poi &b) {
        Type dx=a.x-b.x, dy=a.y-b.y;
        return dx*dx + dy*dy;
    }
    static int ccw(const Poi &a, const Poi &b, const Poi &c) {
        Type d=a.x*b.y + b.x*c.y + c.x*a.y - a.x*c.y - b.x*a.y - c.x*b.y;
        return (d>0)-(d<0);
    }
};

/* should not use
template <typename Type>
struct Line {
    constexpr long double EPS=1e-12L;
    Poi<Type> a, b;
    Line()=default;
    explicit Line(Type m, Type c): a(0, c), b(1, m+c) {}
    explicit Line(const Poi<Type> &a, const Poi<Type> &b): a(a), b(b) {}

    pair<int, Poi<Type>> intersect(const Line &o) const {
        using ldb=long double;
        Poi v1=b-a, v2=o.b-o.a;
        ldb det=v1/v2;

        auto isPoint=[](const Poi<Type> &v) -> bool {
            return v==Poi<Type>(0, 0);
        };
        bool p1=isPoint(v1), p2=isPoint(v2);

        if (p1 && p2) return a==o.a? pair(1, a): pair(0, Poi<Type>());
        if (p1) return fabsl(v2/(a-o.a))<=EPS? pair(1, a)  : pair(0, Poi<Type>());
        if (p2) return fabsl(v1/(o.a-a))<=EPS? pair(1, o.a): pair(0, Poi<Type>());

        if (fabsl(det)<=EPS) {
            ldb cross=(o.a-a)/v1;
            if (fabsl(cross)<=EPS)
                return {2, a};
            return {0, Poi<Type>()};
        }

        ldb num=(o.a-a)/v2;
        ldb t=num/det;
        Poi<Type> p=a+t*v1;
        return {1, p};
    }
};
*/

template <typename Type>
vector<Poi<Type>> convex(vector<Poi<Type>> pts) {
    sort(pts.begin(), pts.end());
    pts.erase(unique(pts.begin(), pts.end()), pts.end());
    if (pts.size()<3) return pts;
    vector<Poi<Type>> H;
    H.reserve(pts.size()*2);

    for (int i=0; i<2; i++) {
        int base=H.size();
        for (auto &p: pts) {
            while (H.size()>base+1 && ccw(H[H.size()-2], H.back(), p)<=0)
                H.pop_back();
            H.push_back(p);
        }
        H.pop_back();
        reverse(pts.begin(), pts.end());
    }
    if (!H.empty() && H.front()==H.back())
        H.pop_back();
    return H;
}

/* TODO:
 * diameter(const vector<Poi> &H)
 * width(const vector<Poi> &H)
 * has many error for Line struct
 */
