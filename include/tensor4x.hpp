#ifndef _TENSOR4X_HPP_
#define _TENSOR4X_HPP_

#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <map>
#include <tuple>
#include <vector>
#include <set>

using namespace std;

typedef tuple<size_t,size_t,size_t,size_t> Index4i_t;

#define cpind(i,j) ((i>j)?(i*(i+1)/2+j):(j*(j+1)/2+i))
#define UNPACK(ind) make_tuple(get<0>(ind),get<1>(ind),get<2>(ind),get<3>(ind))

set<Index4i_t> _symmetrize_index_ (const Index4i_t& ind);

template <class T>
class TENSOR4X
{
    public:
        map<Index4i_t, T> tensor;

        TENSOR4X<T> () : TENSOR4X<T> (0,0,0,0) {}
        TENSOR4X<T> (size_t i) : TENSOR4X<T> (i,i,i,i) {}
        TENSOR4X<T> (size_t i, size_t j, size_t k, size_t l) :
            dim (vector<size_t>({i,j,k,l})), prec (6), symmetry (true) {}
        TENSOR4X<T> (const size_t K, const vector<T>& V);

//  overload parentheses
        inline
        T& operator()(size_t i, size_t j, size_t k, size_t l)
        {
            return tensor[_canonical_order_ (i,j,k,l)];
        }
        inline
        T operator()(size_t i, size_t j, size_t k, size_t l) const
        {
            try
            {
                return tensor.at (_canonical_order_ (i,j,k,l));
            }
            catch (const out_of_range& e)
            {
                return T(0);
            }
        }

//  overload += and -=
        TENSOR4X<T>& operator+=(const TENSOR4X<T>& rhs);
        TENSOR4X<T>& operator-=(const TENSOR4X<T>& rhs);

//  overload scalar multiplication and quotation
        TENSOR4X<T>& operator*=(const double rhs);
        TENSOR4X<T>& operator/=(const double rhs);

//  precision
        inline int _get_precision_ () const {return prec;}
        inline void _set_precision_ (int p) {prec = p;}

//  dim
        inline size_t _get_dim_ (int i) const {return dim[i];}

//  symmetry
        inline bool _get_symmetry_ () const {return symmetry;}
        inline void _set_symmetry_ (bool sym) {symmetry = sym;}

//  write to path
        void _write_to_file_ (const string& path, string fname = "V") const;
//  io
        void _dump_mat_ (const int Switch, const string& info,
            const int B = 0) const;

    private:
        vector<size_t> dim;     // dimension
        int prec;               // cout precision (default: 6)
        bool symmetry;

        void _erase_zeros_ ()
        {
            for (auto it = tensor.cbegin(); it != tensor.cend();)
            {
                if (fabs (it->second) < 1.E-10) tensor.erase(it++);
                else                            ++it;
            }
        };
        void _test_range_ (size_t, size_t, size_t, size_t) const;
        void _test_dim_ (const TENSOR4X<T>& t1, const TENSOR4X<T>& t2);
        Index4i_t _canonical_order_ (size_t, size_t, size_t, size_t);
        Index4i_t _canonical_order_ (size_t, size_t, size_t, size_t) const;
};

//  construct from std::vector<T> 
template <class T>
TENSOR4X<T>::TENSOR4X (const size_t K, const vector<T>& V) :
	TENSOR4X<T> (K)
{
    for (int i=0, ij=0; i < K; i++) for (int j = 0; j <= i; j++, ij++)
    for (int k=0, kl=0; k < K; k++) for (int l = 0; l <= k; l++, kl++)
        if (kl <= ij && fabs (V[cpind(ij,kl)]) > 1.E-10)
            tensor[_canonical_order_ (i,j,k,l)] = V[cpind(ij,kl)];
}

//  overload += and -=
template <class T>
TENSOR4X<T>& TENSOR4X<T>::operator+=(const TENSOR4X<T>& rhs)
{
    _test_dim_ (*this, rhs);
    const map<Index4i_t,T> &rt = rhs.tensor;
    for (auto iter = rt.begin (); iter != rt.end (); ++iter)
        tensor[iter->first] += iter->second;
//  erase zero elements
    _erase_zeros_ ();
    return *this;
}

template <class T>
TENSOR4X<T>& TENSOR4X<T>::operator-=(const TENSOR4X<T>& rhs)
{
    _test_dim_ (*this, rhs);
    const map<Index4i_t,T> &rt = rhs.tensor;
    for (auto iter = rt.begin (); iter != rt.end (); ++iter)
        tensor[iter->first] -= iter->second;
//  erase zero elements
    _erase_zeros_ ();
    return *this;
}

//  overload + and -
template <class T>
TENSOR4X<T> operator+(TENSOR4X<T> lhs, const TENSOR4X<T>& rhs)
{
    lhs += rhs;
    return lhs;
}

template <class T>
TENSOR4X<T> operator-(TENSOR4X<T> lhs, const TENSOR4X<T>& rhs)
{
    lhs -= rhs;
    return lhs;
}

//  overload scalar multiplication and quotation
template <class T>
TENSOR4X<T>& TENSOR4X<T>::operator*=(const double rhs)
{
    for (auto iter = tensor.begin (); iter != tensor.end (); ++iter)
        iter->second *= rhs;
}

template <class T>
TENSOR4X<T>& TENSOR4X<T>::operator/=(const double rhs)
{
    this *= 1. / rhs;
}

template <class T>
TENSOR4X<T> operator*(TENSOR4X<T> lhs, const double rhs)
{
    lhs *= rhs;
    return lhs;
}

template <class T>
TENSOR4X<T> operator*(const double lhs, TENSOR4X<T> rhs)
{
    rhs *= lhs;
    return rhs;
}

template <class T>
TENSOR4X<T> operator/(TENSOR4X<T> lhs, const double rhs)
{
    lhs /= rhs;
    return lhs;
}

template <class T>
TENSOR4X<T> operator/(const double lhs, TENSOR4X<T> rhs)
{
    rhs /= lhs;
    return rhs;
}

//  overload ostream
template <class T>
string to_string_prec (const T val, int prec)
{
    std::ostringstream out;
    out << std::setprecision(prec) << val;
    return out.str();
}

template <class T>
string _to_string_ (const Index4i_t& ind, T val, int prec)
{
    size_t i, j, k, l;
    tie (i, j, k, l) = UNPACK(ind);
    return to_string(i) + ";" + to_string(j) + ";"
        + to_string(k) + ";" + to_string(l) + ";"
        + to_string_prec(val, prec) + "\n";
}

template <class T>
void _symmetrize_strings_ (ostream& os, const Index4i_t& ind, T val, int prec)
{
    set<Index4i_t> ind_set = _symmetrize_index_ (ind);
    set<Index4i_t>::iterator iter;
    for (iter = ind_set.begin (); iter != ind_set.end (); ++iter)
        os << _to_string_ (*iter, val, prec);
}

template <class T>
ostream& operator<<(ostream& os, const TENSOR4X<T>& obj)
{
    os.precision (obj._get_precision_ ());
    for (auto iter = obj.tensor.begin ();
        iter != obj.tensor.end (); ++iter)
    {
        if (obj._get_symmetry_ ())
            os << _to_string_ (iter->first, iter->second, obj._get_precision_ ());
        else
            _symmetrize_strings_ (os, iter->first, iter->second,
                obj._get_precision_ ());
    }
    return os;
}

/*  using the following eigth-fold symmetry
*       (ij|kl) = (ji|kl) = (ij|lk) = (ji|lk)
*     = (kl|ij) = (kl|ji) = (lk|ij) = (lk|ji)
*   to put input index set in to the canonical order:
*       j <= i, l <= k and cpind(k,l) <= cpind(i,j)
*/
template <class T>
Index4i_t TENSOR4X<T>::_canonical_order_ (
    size_t i, size_t j, size_t k, size_t l) const
{
    _test_range_ (i,j,k,l);
    if (j > i)  swap (i, j);
    if (l > k)  swap (k, l);
    size_t kl = cpind(k,l), ij = cpind(i,j);
    if (kl > ij)
        tie (i,j,k,l) = make_tuple (k,l,i,j);
    return make_tuple (i,j,k,l);
}

template <class T>
Index4i_t TENSOR4X<T>::_canonical_order_ (
    size_t i, size_t j, size_t k, size_t l)
{
    _test_range_ (i,j,k,l);
    if (j > i)  swap (i, j);
    if (l > k)  swap (k, l);
    size_t kl = cpind(k,l), ij = cpind(i,j);
    if (kl > ij)
        tie (i,j,k,l) = make_tuple (k,l,i,j);
    return make_tuple (i,j,k,l);
}


//  test_range
template <class T>
void TENSOR4X<T>::_test_range_ (size_t i, size_t j, size_t k, size_t l) const
{
    try
    {
        if(i >= dim[0] || j >= dim[1] || k >= dim[2] || l >= dim[3])
        {
            string err = "[ERROR] index requested ("
                + to_string(i) + "," + to_string(j) + ","
                + to_string(k) + "," + to_string(l)
                + ") is out of range in TENSOR4X.";
            throw out_of_range (err);
        }
    }
    catch (const out_of_range& e)
    {
        cout << e.what () << endl;
        exit (1);
    }
}

//  test dimension
template <class T>
void TENSOR4X<T>::_test_dim_ (const TENSOR4X<T>& t1, const TENSOR4X<T>& t2)
{
    try
    {
        if (t1._get_dim_ (0) != t2._get_dim_ (0) ||
            t1._get_dim_ (1) != t2._get_dim_ (1) ||
            t1._get_dim_ (2) != t2._get_dim_ (2) ||
            t1._get_dim_ (3) != t2._get_dim_ (3))
            throw range_error ("[ERROR] dimension must match"
                " for two TENSOR4X-type variables to be operated.");
    }
    catch (const range_error& e)
    {
        cout << e.what () << endl;
        exit (1);
    }

}

//  write to path
template <class T>
void TENSOR4X<T>::_write_to_file_ (const string& path, string fname) const
{
    string fnameV (path + "/" + fname);
    FILE * pV = fopen (fnameV.c_str (), "w+");
    for (auto iter = tensor.cbegin (); iter != tensor.cend (); ++iter)
        if (fabs (iter->second) > 1.E-10)
        {
            size_t i,j,k,l; tie (i,j,k,l) = UNPACK(iter->first);
            fprintf (pV, "%d;%d;%d;%d;%18.16f\n", i,k,j,l,iter->second);
        }
    fclose (pV);
}

template <class T>
void TENSOR4X<T>::_dump_mat_ (const int Switch, const string& info,
    const int B) const
{
    cout << info << endl;
    for (int mu = 0; mu < dim[0]; mu++)
    {
        for (int nu = 0; nu < dim[0]; nu++)
            if (Switch == 1)
                printf ("%8.5f ", (*this)(mu,mu,nu,nu));
            else if (Switch == 2)
                printf ("%8.5f ", (*this)(mu,nu,mu,nu));
            else if (Switch == 3)
            {
                if (mu < nu)
                    printf ("%8.5f ", (*this)(mu,nu,nu,nu));
                else
                    printf ("%8.5f ", (*this)(nu,mu,mu,mu));
            }
            else if (Switch == 4)   // BBCD
                printf ("%8.5f ", (*this)(B,B,mu,nu));
            else if (Switch == 5)   // BCBD
                printf ("%8.5f ", (*this)(B,mu,B,nu));
        cout << endl;
    }
    cout << endl;
}

#endif
