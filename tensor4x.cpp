#include "include/tensor4x.hpp"

set<Index4i_t> _symmetrize_index_ (const Index4i_t& ind)
{
	size_t i, j, k, l, cpij, cpkl;
    tie (i, j, k, l) = UNPACK(ind); cpij = cpind(i,j);  cpkl = cpind(k,l);
    bool ij = (i == j), kl = (k == l), ijkl = (cpij == cpkl);
    set<Index4i_t> ind_set; ind_set.insert (ind);
    set<Index4i_t>::iterator iter;
	if (!ij)
        for (iter = ind_set.begin (); iter != ind_set.end (); ++iter)
        {
            tie (i, j, k, l) = UNPACK(*iter);
            swap (i, j);
            ind_set.insert (make_tuple (i, j, k, l));
        }
    if (!kl)
        for (iter = ind_set.begin (); iter != ind_set.end (); ++iter)
        {
            tie (i, j, k, l) = UNPACK(*iter);
            swap (k, l);
            ind_set.insert (make_tuple (i, j, k, l));
        }
    if (!ijkl)
        for (iter = ind_set.begin (); iter != ind_set.end (); ++iter)
        {
            tie (i, j, k, l) = UNPACK(*iter);
            swap (i, k);    swap (j, l);
            ind_set.insert (make_tuple (i, j, k, l));
        }
	return ind_set;
}

