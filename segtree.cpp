#include <bits/stdc++.h>
using namespace std;

template<typename Type, typename M>
concept Mergable=requires(M merge, Type a, Type b) {
    { merge(a, b) } -> same_as<Type>;
    // Type() must be an identity for merge()
};

template<typename Type, typename M>
requires Mergable<Type, M>
class SegTree {
private:
    size_t n;
    vector<Type> tree;
    M merge;

    void set(size_t idx, Type val) {
        for (tree[idx+=n]=val; idx>1; idx>>=1)
            tree[idx/2]=merge(tree[idx], tree[idx^1]);
    }

public:
    SegTree(size_t n): n(n) {
        tree=vector<Type>(n*2);
    }
    SegTree(const vector<Type> &arr): n(arr.size()) {
        tree=vector<Type>(n*2);
        for (size_t i=0; i<n; i++)
            tree[i+n]=arr[i];
        for (size_t i=n-1; i>0; i--)
            tree[i]=merge(tree[i*2], tree[i*2+1]);
    }

    struct Proxy {
        SegTree &tree;
        size_t idx;

        Proxy &operator=(Type val) {
            tree.set(idx, val);
            return *this;
        }

        operator Type() const {
            return tree.tree[idx+tree.n];
        }
    };
    Proxy operator[](size_t idx) {
        return Proxy(*this, idx);
    }

    Type query(size_t l, size_t r) {
        Type k=Type();
        for (l+=n, r+=n; l<=r; l>>=1, r>>=1) {
            if (l&1)
                k=merge(k, tree[l++]);
            if (~r&1)
                k=merge(k, tree[r--]);
        }
        return k;
    }
};

class Add {
public:
    int operator()(int x, int y) {
        return x+y;
    }
};

int main() {
    cin.tie(0)->sync_with_stdio(0);
    SegTree<int, Add> tree(vector<int>(10, 1));

    cout << tree.query(0, 3) << '\n';
    tree[5]=2, tree[6]=7;
    cout << tree.query(4, 7) << '\n';
}