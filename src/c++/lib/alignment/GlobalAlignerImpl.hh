
#include <algorithm>

#ifdef ALN_DEBUG
#include <iostream>
#endif


template <typename ScoreType>
template <typename SymIter>
void
GlobalAligner<ScoreType>::
align(
    SymIter begin1, SymIter end1,
    SymIter begin2, SymIter end2,
    AlignmentResult& result)
{
#ifdef ALN_DEBUG
    std::ostream& log_os(std::cerr);
#endif

    result.apath.clear();

    const size_t size1(std::distance(begin1, end1));
    const size_t size2(std::distance(begin2, end2));
    _scoreMat.resize(size1+1, size2+1);
    _whichMat.resize(size2+1, size2+1);

    for (unsigned i(0); i<=size1; i++)
    {
        E_[i][0]=-10000;
        F_[i][0]=-10000;
        _scoreMat.match[i][0]=-10000;
    }

    for (unsigned j(0); j<=y_size; j++)
    {
        F_[0][j]=0;
        E_[0][j]=-10000;
        _scoreMat.match[0][j]=0;
    }
    for (unsigned i(1); i<=x_size; i++)
    {
#ifdef ALN_DEBUG
        log_os << i << ": ";
#endif
        for (unsigned j(1); j<=y_size; j++)
        {
            // update _scoreMat.match
            max3(_scoreMat.match[i][j],T_scoreMat.match[i][j],
                 _scoreMat.match[i-1][j-1],E_[i-1][j-1] - w_open_, F_[i-1][j-1] - w_open_);

            if (x[i-1]==y[j-1])
            {
                _scoreMat.match[i][j]+=w_match_;
            }
            else
            {
                _scoreMat.match[i][j]-=w_mismatch_;
            }

            // update E_
            max3(E_[i][j],TE_[i][j],
                 _scoreMat.match[i][j-1]-w_open_,
                 E_[i][j-1],
                 F_[i][j-1]-w_open_);

            E_[i][j]-=w_extend_;

            // update F_
            max3(F_[i][j],TF_[i][j],
                 _scoreMat.match[i-1][j]-w_open_,
                 E_[i-1][j]-w_open_,
                 F_[i-1][j]);

            F_[i][j]-=w_extend_;
#ifdef ALN_DEBUG
            printf(" %3d:%3d:%3d/%1d%1d%1d",_scoreMat.match[i][j],E_[i][j],F_[i][j],T_scoreMat.match[i][j],TE_[i][j],TF_[i][j]);
#endif
        }
#ifdef ALN_DEBUG
        log_os << "\n";
#endif
    }

    ScoreType max(0),thisMax;
    uint8_t whichMatrix(0),thisMatrix;
    int ii(0),jj(0);

    for (unsigned j(0); j<=y_size; j++)
    {
        if (_scoreMat.match[x_size][j] > F_[x_size][j]) {
            thisMax=_scoreMat.match[x_size][j];
            thisMatrix=0;
        } else {
            thisMax=F_[x_size][j];
            thisMatrix=2;
        }

        if ((j==0) || (thisMax>max))
        {
            max=thisMax;
            ii=x_size;
            jj=j;
            whichMatrix=thisMatrix;
        }
    }


    result.myScore = max;
#ifdef ALN_DEBUG
    cout << "start cell = " << ii << " " << jj << " " << whichMatrix << endl;
#endif
    result.x_start_ = ii;
    result.y_start_ = jj;

    DPMatrixN* pT_ =&TE_;

    while ((ii>0)&&(jj>0))
    {
        if (whichMatrix==0)
        {
            pT_=&T_scoreMat.match;
        }
        else if (whichMatrix==1)
        {
            pT_=&TE_;
        }
        else
        {
            //assert(whichMatrix==2);
            pT_=&TF_;
        }

#ifdef ALN_DEBUG
        log_os << ii << " " << jj << " " << whichMatrix << " "
               << " "<< _scoreMat.match[ii][jj]<< ":"<< E_[ii][jj] << ":" << F_[ii][jj]
               << "/"
               << " "<< T_scoreMat.match[ii][jj] << ":"<< TE_[ii][jj] << ":" << TF_[ii][jj]
               << " - " << (*pT_)[ii][jj] << "\n"
#endif
        const uint8_t nextMatrix=(*pT_)[ii][jj];
        if (whichMatrix==0)
        {
            xt_+=x[ii-1];
            yt_+=y[jj-1];
            at_+=((x[ii-1]==y[jj-1])?'|':'.');
            ii--;
            jj--;
            // whichMatrix=0;
        } // ~if
        else if (whichMatrix==1)
        {
            xt_+='-';
            yt_+=y[jj-1];
            at_+=' ';
            jj--;
            //    whichMatrix=1;
        } // ~else if
        else
        {
            //assert(whichMatrix==2);
            xt_+=x[ii-1];
            yt_+='-';
            at_+=' ';
            ii--;
            //      whichMatrix=2;
        } // ~else
        whichMatrix=nextMatrix;
    } // ~while
#ifdef ALN_DEBUG
    log_os << "end cell = " << ii << " " << jj << "\n";
#endif
    result.x_end_ = ii;
    result.y_end_ = jj;

    reverse(xt_.begin(),xt_.end());
    reverse(at_.begin(),at_.end());
    reverse(yt_.begin(),yt_.end());
}

