#include "mmC.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


SEXP DPre(SEXP yvec, SEXP Xmat, SEXP Zmat,  SEXP personsvec, SEXP niterInt,
        SEXP nburnInt, SEXP nthinInt, SEXP shapealph, SEXP ratebeta)
{
BEGIN_RCPP
    // Time run
    //wall_clock timer;
    //timer.tic();
    // Wrap SEXP objects in Rcpp Vector and Matrix objects.  Not copies,
    // but wrappers.
    // Data: observations, fe, re, mm weights, car weights, session groups,
    // client id's
    // Counters (number of); cbt cases, mcmc iter, burnin.
    NumericVector yr(yvec);
    NumericMatrix Xr(Xmat);
    NumericMatrix Zr(Zmat);
    IntegerVector pr(personsvec);
    int niter = as<int>(niterInt);
    int nburn = as<int>(nburnInt);
    int nthin = as<int>(nthinInt);
    int nkeep = (niter -  nburn)/nthin;
    double ac = as<double>(shapealph);
    double bc = as<double>(ratebeta);

    // Extract key row and column dimensions we will need to sample parameters
    int nc =  Xr.nrow(), nf = Xr.ncol(), nr = Zr.ncol();

    // Compute np, number of unique clients
    // algorithm assumes cases clustered by unique client
    IntegerVector dpr = diff(pr);
    icolvec diffpersons(dpr.begin(), nc - 1, false);
    int np = accu(diffpersons) + 1;
    /* int np = accumulate(dpr.begin,dpr.end,1); */


    // Create armadillo structures we will use for our computations
    // We set copy_aux_mem = false from the constructor, meaning we re-use
    // the memory to in which our object is stored
    // Xr.begin() is a pointer, presumably of type double*
    mat xmat(Xr.begin(), nc, nf, false);
    mat zmat(Zr.begin(), nc, nr, false);
    colvec y(yr.begin(), nc, false);
    icolvec persons(pr.begin(), nc, false);
    field<mat> zsplit(np,1);
    field<colvec> ytilsplit(np,1);

    // Set random number generator state
    RNGScope scope; /* Rcpp */
    srand ( time(NULL) ); /* arma */

    // Armadillo structures to capture results
    // Will be wrapped into a list element on return
    colvec Deviance(nkeep); Deviance.zeros();
    colvec Alpha(nkeep);
    colvec Conc(nkeep);
    mat Devmarg(nkeep,nc);
    mat Resid(nkeep,nc);
    mat Beta(nkeep,nf);
    mat B(nkeep,nr*np);
    rowvec brow(np*nr); brow.zeros();
    mat Taub(nkeep,nr);
    colvec Taue(nkeep);
    icolvec numM(nkeep);
    imat S(nkeep,np);
    field<icolvec> Num(nkeep,1);
    field<ucolvec> bigS;
    field<mat> optPartition(3,1); /* Hold empirical probability of co-clustering matrix, index of L-sq clustering scores objects, and cluster identifiers per iteration */

    // Initialize parameter values
    /* precision parameters */
    double taubeta, taualph, taue;
    taubeta = taualph = taue = 1;
    rowvec taub(nr); taub.ones();
    /* cluster capture variables */
    IntegerVector ss = seq_len(np);
    icolvec s(ss.begin(), np, false);
    s -= 1; s = shuffle(s);
    int M = np; /* each person initialized with their own cluster */
    icolvec num(M); num.ones();
    mat    bstarmat = randn<mat>(M,nr)*sqrt(1/taub(0));
    mat    bmat(np,nr); bmat.zeros();
    colvec zb(nc); zb.zeros();
    double conc = 1; /* DP concentration parameter */
    /* remaining parameters */
    colvec beta = randn<colvec>(nf)*sqrt(1/taubeta);
    colvec resid(nc); resid.zeros(); /* regression residual */
    double alpha = rnorm( 1, 0, sqrt(1/taualph) )[0];
    mat pbmat(nr,nr); pbmat.zeros();
    /* fit assessment - related measures */
    double deviance = 0;  ucolvec ordscore; ordscore.zeros();
    mat phat(np,np); phat.zeros(); /* empirical co-clustering probability matrix */
    rowvec devmarg(nc); colvec devres(4); devres.zeros();
    colvec devres8(4); devres8.zeros();
    rowvec logcpo(nc); logcpo.zeros(); double lpml;
    /* capture samples to return */
    int oo = 0, kk;

    // set hyperparameter values
    double a2, a3, a6, b2, b3, b6;
    a2 = b2 = 1.0; /* taub */
    a3 = b3 = 1; /*taue */
    a6 = ac; b6 = bc; /* conc */
    
    // Conduct posterior samples and store results
    int k,l;
    for(k = 0; k < niter; k++)
    {
        //if( (k % 1000) == 0 ) Rcout << "Interation: " << k << endl;
        clusterstep(xmat, zmat, zsplit, ytilsplit, pbmat, y, beta, alpha,
                taue, taub, persons, bstarmat, s, num, M, conc,
                np, nr);
        concstep(conc, M, np, a6, b6);
        bstarstep(zsplit, ytilsplit, taue, pbmat, s, num, bstarmat,
                zb, bmat, np, nr);
        betaDPstep(xmat, y, beta, alpha, taue, zb, nf);
        alphaDPstep(xmat, beta, zb, y, resid, alpha, taue, nc);
        tauDPstep(bstarmat, resid, taub, taue, a2, a3, b2, b3, M, nc);
        if(k >= nburn)
        {
            deviance =  dev(resid, taue); /* scalar double deviance */
            dmarg(resid, taue, devmarg); /* 1 x nc vector of densities*/
            kk = k - nburn;
            if( kk == ((oo+1)*nthin - 1) )
            {
                Deviance(oo) = deviance;
                Devmarg.row(oo)  = devmarg;
                Beta.row(oo) = trans(beta);
                Alpha(oo) = alpha;
                Conc(oo) = conc;
                for( l = 0; l < nr; l++)
                {
                brow.cols( np*l,(np*(l+1)-1) ) = trans( bmat.col(l) );
                }
                B.row(oo) = brow;
                Taue(oo) = taue;
                Taub.row(oo) = taub;
                Resid.row(oo) = trans(resid);
                numM(oo) = M;
                S.row(oo) = trans(s);
                Num(oo,0) = num;
                oo += 1; /* increment sample return counter */
            } /* end conditional statement on returning sample */
        } /* end if k > burnin, record results for return */
    } /* end MCMC loop over k */

    // compute least squares cluster
    lsqcluster(S, Num, ordscore, phat, bigS);
    optPartition(0,0) = phat;
    optPartition(1,0) = conv_to<mat>::from(ordscore);
    optPartition(2,0) = conv_to<mat>::from(S);
    // DIC
    dic3comp(Deviance, Devmarg, devres); /* devres = c(dic,dbar,dhat,pd) */
    dic8dpcomp(Deviance, Alpha, Beta, B, Taue, xmat, zmat,
            y, persons, devres8, np);
    cpo(Devmarg, logcpo, lpml);

    // Return results
    return Rcpp::List::create(Rcpp::Named("Deviance") = Deviance,
                                  Rcpp::Named("Devmarg") = Devmarg,
                                  Rcpp::Named("devres") = devres,
                                  Rcpp::Named("devresp") = devres8,
                                  Rcpp::Named("logcpo") = logcpo,
                                  Rcpp::Named("lpml") = lpml,
				  Rcpp::Named("Beta") = Beta,
				  Rcpp::Named("Alpha") = Alpha,
                                  Rcpp::Named("C") = Conc,
                                  Rcpp::Named("B") = B,
                                  Rcpp::Named("Taue") = Taue,
                                  Rcpp::Named("Taub") = Taub,
                                  Rcpp::Named("Residuals") = Resid,
                                  Rcpp::Named("M") = numM,
                                  //Rcpp::Named("S") = S,
                                  Rcpp::Named("Num") = Num,
                                  Rcpp::Named("optPartition") = optPartition,
                                  Rcpp::Named("bigSmin") = bigS
				  );

    // Print CPU runtime
     //double n_secs = timer.toc();
     //cout << "took " << n_secs << " seconds" << endl;
END_RCPP
} /* end MCMC function returning SEXP */


    // static function to sample random effects
    // updates bmat (npxnr) and zb[{i:person(i)==j}] = {zj * bmat[j,]}
    SEXP clusterstep(const mat& xmat, const mat& zmat, field<mat>& zsplit,
            field<colvec>& ytilsplit, mat& pbmat, const colvec& y, 
            const colvec& beta, double alpha, double taue, const rowvec& taub,
            icolvec& persons, mat& bstarmat,
            icolvec& s, icolvec& num, int& M, double& conc,
            int np, int nr)
    {
        BEGIN_RCPP
        // compute number of cases for each person and store results in
        // 'positions' (np x 1)
        icolvec positions(np); positions.zeros();
        icolvec::iterator aa = persons.begin();
        icolvec::iterator bb = persons.end();
        for(icolvec::iterator i = aa; i != bb ; i++)
        {
           positions(*i-1) += 1;
        }

        // sample cluster assignments, s(1), ..., s(np)
        // fixing all other parameters, including bstarmat
        mat zj, xj, yj, phib(nr,nr);
        colvec cj, cstarj, ytildej, bstarj(nr), eb(nr), hb(nr);
        double logq0, q0, sweights;
        int startrow = 0; 
        int j, l, m, nj, endrow;
        for(j = 0; j < np; j++)
        {
            // store by-person data cuts for later sampling of clust locs
            /* extract by-person matrix views from data*/
            nj = positions(j);
            endrow = startrow + nj - 1;
            zj = zmat.rows(startrow,endrow);
            xj = xmat.rows(startrow,endrow);
            yj = y.rows(startrow,endrow);
            zsplit(j,0) = zj; /* np rows */

            /* subtracting cj eliminates need to t-fer xj to sample bstar */
            cj = alpha + xj*beta;
            ytildej = yj - cj;
            ytilsplit(j,0) = ytildej;

            // sample cluster assignment indicators, s(1),...,s(np)
            if(num(s(j)) == 1) /* remove singleton cluster */
            {
                bstarmat.shed_row(s(j));
                num.shed_row(s(j));
                M -= 1; /* decrement cluster count */

                //decrement cluster tracking values by 1 for tossed cluster
                uvec ord = sort_index(s); /* of size np, same as s */
                int k; int thresh = np;
                for (k = 0; k < np; k++)
                {
                    if( s(ord(k)) > s(j) )
                    {
                        thresh = k;
                        break;
                    }
                } /* end loop k to find thresh*/
                // decrement s > s(j)
                for (m = thresh; m < np; m++)
                {
                    s(ord(m)) -= 1;
                }
            } /* end cluster accounting adjustment for singleton cluster */
            else
            {
                num(s(j)) -= 1;
            } /* decrement non-singleton cluster count by one */

            // construct normalization constant, q0, to sample s(j)
            // compute posterior mean and variance to compute density for bj
            // build loqq0 and exponentiate
            pbmat.eye(nr,nr);
            NumericVector zro(1);
            logq0 = loglike(ytildej,taue);
            int i;
            for(i = 0; i < nr; i++)
            {
                pbmat.diag()(i) = taub(i);
                logq0 += dnorm(zro, 0.0, sqrt(1/taub(i)), true)[0]; /* true = log */
            }
            eb = taue*trans(zj)*ytildej;
            phib = taue*trans(zj)*zj + pbmat; 
	    hb = solve(phib,eb);
            logq0 -= logdens(hb,phib); /* dmvn(hb,phib^-1) */
            q0 = exp(logq0);

            // construct posterior sampling weights to sample s(j)
            colvec weights(M+1); weights.zeros();
            for (l = 0; l < M; l++) /* cycle through all clusters for s(j) */
            {
                s(j) = l; /* will compute likelihoods for every cluster */
                bstarj = trans( bstarmat.row(s(j)) );
                cstarj = cj + zj*bstarj;
                ytildej = yj - cstarj;

                weights(l) = exp(loglike(ytildej,taue));
                weights(l) *= num(s(j))/(double(np) - 1 + conc);
            }
            weights(M) = conc/(double(np)-1+conc)*q0;

            // normalize weights
            sweights = sum(weights);
            if(sweights == 0)
            {
                weights.ones(); weights *= 1/(double(M)+1);
            }
            else
            {
                weights = weights / sweights;
            }

            // conduct discrete posterior draw for s(j)
            unsigned long MplusOne = M + 1;
            s(j) = rdrawone(weights, MplusOne);

            // if new cluster chosen, generate new location
            if(s(j) == M)
            {
                bstarj.zeros();
                //rmvnqr(zj, Ub, yj, cj, bstarj, nj, nr, taue);
                rmvnchol(zj, pbmat, yj, cj, bstarj, nr, taue);
                bstarmat.insert_rows(M,1);
                num.insert_rows(M,1);
                bstarmat.row(M) = trans(bstarj);
                num(M) = 1;
                M = MplusOne;
            }
            else
            {
                num(s(j)) += 1;
            }
            
            //re-set start positions
            startrow += nj;
        } /* end loop j for cluster assignment */
        END_RCPP
    } /* end function bstep for cluster assignments, s, and computing zb */

    SEXP concstep(double& conc, int M, int np,  double a6, double b6)
    {
        BEGIN_RCPP
        // sample DP concentration parameter, c
        // c|s independent of sample, y
        // see page 585 of Escobar and West (1995)
        double eta = rbeta( 1,conc+1,double(np) )[0];
        double drawP;
        drawP = (a6+double(M)-1)/( a6+double(M)-1 + double(np)*(b6-log(eta)) );
        int mixcomp = rbinom(1,1,drawP)[0];
        if( mixcomp == 1 )
        {
            conc = rgamma( 1, a6 + M, 1/(b6-log(eta)) )[0];
        }
        else
        {
            conc = rgamma( 1, a6 + M - 1, 1/(b6-log(eta)) )[0];
        }
        END_RCPP
    } /* end function to sample concentration parameter, conc */
    
    SEXP bstarstep(const field<mat>& zsplit, const field<colvec>& ytilsplit,
            double taue, const mat& pbmat,
            const icolvec& s, const icolvec& num, mat& bstarmat,
            mat& zb, mat& bmat, int np, int nr)
    {
        BEGIN_RCPP
        // sample bmat, independently, by person
        int M = bstarmat.n_rows;
        int j, k, l, m, numj;
        colvec bstarj(nr);
        
        // collect subsets of cluster indices, i: s(i) = j, j = 1, ..., M

        uvec ord = sort_index(s); /* of length np > = M*/
        uvec thresh(M); thresh(0) = 0;
        int count = 1;
        //find cluster index start (change) points
        for (m = 0; m < (np-1); m++)
        {
            if( s(ord(m+1)) > s(ord(m)) )
            {
                thresh(count) = m+1;
                count += 1;
            }
        }
        // sample bstarmat for each cluster, j
        int rowchoose;
        for(j=0; j < M; j++) /* j is now the cluster label */
        {
            numj = num(j);
            field<mat> zjwedge(numj,1);
            field<colvec> yjtilwedge(numj,1);
            // collect data observations - z,ytilde - for each cluster
            for(k = 0; k < numj; k++)
            {
                /* extract by-person matrix views from data*/
                rowchoose = ord(thresh(j) + k);
                zjwedge.row(k) = zsplit.row(rowchoose); /* zsplit is a field<mat> containing repeated obs for person = rowchoose */
                yjtilwedge.row(k) = ytilsplit.row(rowchoose);
            }
            /* sample posterior of nj block of bmat*/
            bstarj.zeros();
            //rmvnqrclust(zjwedge, Ub, yjtilwedge, bstarj, nr, taue);
            rmvncholclust(zjwedge, pbmat, yjtilwedge, bstarj, taue);
            bstarmat.row(j) = trans(bstarj);

        } /* end loop j for sampling bstarj */

        // construct bmat for function return
        // compute zbl (lth piece) of zb = [z1*b1vec,...,znp*bnpvec]
         int startrow = 0; int endrow, nl; rowvec blr(nr);
         for(l = 0; l < np; l++) /* now we're iterating by person */
            {
                /* bmat */
                blr = bstarmat.row(s(l));
                bmat.row(l) = blr;
                /* zb */
                nl = zsplit(l,0).n_rows;
                endrow = startrow + nl - 1;
                zb.rows(startrow,endrow) = zsplit(l,0)*trans(blr);
                startrow += nl;
            }

        END_RCPP
    } /* end function bstep for sampling bmat and zb */

    SEXP betaDPstep(const mat& xmat, const colvec& y,
                colvec& beta, double alpha, double taue,
                const colvec& zb, int nf)
    {
        BEGIN_RCPP
        // Sample posterior of beta (nf x 1) from posterior Gaussian
        colvec c = alpha + zb; /* nc x 1 */
        rmvnlschol(xmat, y, c, beta, nf, taue);
        END_RCPP
    } /* end function to sample nf x 1 fixed effects, beta */


    SEXP alphaDPstep(const mat& xmat, const colvec& beta,
                const colvec& zb, const colvec& y,
                colvec& resid, double& alpha, double taue,
                long nc)
    {
        BEGIN_RCPP
        colvec c = xmat*beta + zb;
        resid = y - c;
        double ea = taue*sum(resid);
        double phia = taue*double(nc);
        double ha = ea*(1/phia);
        alpha = rnorm( 1, ha, sqrt(1/phia) )[0];
        resid -= alpha;
        END_RCPP
    } /* end function to sample intercept, alpha */

    SEXP tauDPstep(const mat& bstarmat, const colvec& resid, 
            rowvec& taub, double& taue, double a2, double a3, 
            double b2, double b3, int M, int nc)
    {
        BEGIN_RCPP
        double a, b;
        /* taub */
        int nr = bstarmat.n_cols;
        a = 0.5*double(M) + a2;
        int i;
        for(i = 0; i < nr; i++)
        {
            b =  0.5*dot( bstarmat.col(i), bstarmat.col(i) ) + b2;
            taub(i) = rgamma(1, a, (1/b))[0];
        }
        /* taue */
        a = 0.5*double(nc) + a3;
        b = 0.5*dot( resid, resid ) + b3;
        taue = rgamma(1, a, (1/b))[0];
        END_RCPP
    } /* end function taustep to sample precision parameters */


    double loglike(const colvec& resid, double taue)
    {
        // NumericVector resid.
        // devvec size equal to that of resid.
        NumericVector r = wrap(resid);
        NumericVector devvec  = dnorm( r, 0.0, sqrt(1/taue), true ); /* logD */
        double loglike = accumulate(devvec.begin(),devvec.end(),0.0);
        return(loglike);
    }
    

    double logdens(const colvec& hb, const mat& phib)
    {
        // NumericVector r(resid);
        // devvec is an nc x 1 vector
        double val, sign;
        log_det(val,sign,phib);
        double p = hb.n_rows;
        double c = log(2*M_PI)*(-0.5*p);
 	double logdens = c + 0.5*val - 0.5*as_scalar( trans(hb)*phib*hb );
        return(logdens);
    }


    SEXP dmarg(const colvec& resid, double taue, rowvec& devmarg)
    {
        BEGIN_RCPP
        // NumericVector r(resid);
        // devvec is an nc x 1 vector
        int nc = resid.n_elem;
        NumericVector r = wrap(resid); 
        NumericVector devvec  = dnorm( r, 0.0, sqrt(1/taue), false ); /* f(y|theta)*/
        rowvec devcond(devvec.begin(), nc, false);  devmarg = devcond;
        END_RCPP
    }

    SEXP cpo(const mat& Devmarg, rowvec& logcpo, double& lpml)
    {
        BEGIN_RCPP
        mat invDmarg = pow(Devmarg,-1); /* invert each member of Devmarg */
        logcpo = mean(invDmarg,0); /* 1 x nc */
        logcpo = pow(logcpo,-1); /* invert again to achieve f(y_i|y_-i) */
        logcpo = log(logcpo); /* return vector of log(cpo_i) */
        lpml = sum(logcpo); /* sum each elem of log(cpo_i) */
        END_RCPP
    }
    
    SEXP dic3comp(const colvec& Deviance, const mat& Devmarg, colvec& devres)
    {
        BEGIN_RCPP
        double dbar = mean(Deviance);
        rowvec devmean = mean(Devmarg,0);
        rowvec dens = log(devmean); /* 1 x nc */
        double dhat = sum( dens ); dhat *= -2;
        double dic = 2*dbar - dhat;
        double pd = dic - dbar;
        devres(0) = dic; devres(1) = dbar; devres(2) = dhat; devres(3) = pd;
        END_RCPP
    }

    SEXP dic8dpcomp(const colvec& Deviance, const colvec& Alpha,
            const mat& Beta, const mat& B, const colvec& Taue,
            const mat& xmat, const mat& zmat,
            const colvec& y, icolvec& persons, colvec& devres, int np)
    {
        BEGIN_RCPP
        // compute MAP parameter estimates
        double dbar = mean(Deviance);
        double alpha = mean(Alpha);
        double taue = mean(Taue);
        rowvec beta = mean(Beta,0);
        // compute zbhat
        int nc = xmat.n_rows;
        int nkeep = B.n_rows;
        colvec zb(nc); zb.zeros();
        colvec residp = y - alpha - xmat*trans(beta);
        colvec resid(nc); rowvec dhatrow;
        mat dhatmat(nkeep, nc); dhatmat.zeros();
        int i;
        for( i = 0; i < nkeep; i++ )
        {
            zbcomp(B.row(i), zmat, persons, zb, np);
            resid = residp  - zb;
            dmarg(resid,taue,dhatrow);
            dhatmat.row(i) = dhatrow;
        }
        rowvec dhatvec = mean(dhatmat,0);
        rowvec dens = log(dhatvec);
        double dhat = sum(dens); dhat *= -2;
        double dic = 2*dbar - dhat;
        double pd = dic - dbar;
        devres(0) = dic; devres(1) = dbar; devres(2) = dhat; devres(3) = pd;
        END_RCPP
    }

    unsigned long rdrawone(const colvec& pr, unsigned long k)
    {
        // extract sort index of pr in descending order
        uvec pOrderIndex = sort_index(pr,1);
        // draw uniform random variate we will compare
        // to cdf composed from categories in descending
        // order
        double drawP = runif(1,0,1)[0];
        double pSum = 0.0;
        unsigned long i, x;
        for(i = 0; i < k; i++)
        {
          x = pOrderIndex(i);
          pSum += pr(x);
          if(pSum > drawP)
          {
            return x;
          } /* end conditional test on category i */
        } /* end loop over k categories */
        return pOrderIndex(k-1);
     } /* end function rdrawone */

    
    SEXP rmvnqrclust(const field<mat>& zsplit, const mat& Umat,
            const field<colvec>& ytilsplit, colvec& b, double taue)
    {
        BEGIN_RCPP
        // build xbig and ytildebig
        int p = zsplit(0,0).n_rows;
        double errinv = sqrt(taue);
        int chunks = zsplit.n_rows;
        int i, j, n;
        uvec nchunk(chunks);
        // determine total row size of Xbig
        // by summing over row counts of constituent members
        for(i = 0; i < chunks; i++)
        {
           nchunk(i) += zsplit(i,0).n_rows;
        }
        n = sum(nchunk);
        mat xbig(n+p,p); xbig.zeros();
        colvec ytildebig(n+p); ytildebig.zeros();
        colvec noisevec = randn<colvec>(p);
        // fill xbig, ytildebig
        int startrow = 0;
        int endrow;
        for(j = 0; j < chunks; j++)
        {
            endrow = startrow + nchunk(j) - 1;
            xbig.rows(startrow,endrow) = errinv*zsplit(j,0);
            ytildebig.rows(startrow,endrow) = errinv*ytilsplit(j,0);
            startrow += nchunk(j);
        }
        xbig.rows(n,n+p-1) = Umat;
        // conduct QR backsolving for bcvec
        mat Q, R;
        qr(Q,R,xbig);
        // Translate from Matlab method where
        // for mxn, X where m > n, Q is m x m and R is m x n
        // The m-n cols of Q are basis for OC of X and m-n rows of R are O
        Q = Q.cols(0,p-1);
	R = R.rows(0,p-1);
        colvec q = trans(Q)*ytildebig;
        q += noisevec;
        // return b
        b = solve(R,q);
        END_RCPP
    } /* end function rmvnqrclust for sampling bstarmat */

    SEXP rmvncholclust(const field<mat>& zwedge, const mat& Pmat,
            const field<colvec>& ytilwedge, colvec& b, double taue)
    {
        BEGIN_RCPP
        int n = zwedge.n_rows; int p = zwedge(0,0).n_cols;
        mat ztz(p,p); ztz.zeros(); colvec zty(p); zty.zeros();
        mat phi(p,p);
        int i;
        for(i = 0; i < n; i++)
        {
           ztz += trans(zwedge(i,0))*zwedge(i,0);
           zty += trans(zwedge(i,0))*ytilwedge(i,0);
        }
        phi 		= symmatl(taue*ztz + Pmat); 
	colvec h 	= taue*solve(phi,zty);
        colvec noise 	= randn<colvec>(p);
	b 		= solve(trimatu(chol(phi)),noise) + h;
        END_RCPP
    }

    SEXP lsqcluster(const imat& S, const field<icolvec>& Num,
            ucolvec& ordscore, mat& phat, field<ucolvec>& bigS)
    {
        BEGIN_RCPP
	bigS.set_size(1,1);
        int n = S.n_rows;  int np = S.n_cols;
        mat delta; delta.eye(np,np);
        mat diffsq(np,np); diffsq.zeros();
        irowvec s(np); colvec score(n);
	ucolvec ord(np); ord.zeros();
	unsigned long minscoreloc;
        int i, j, k, l;
        // create phat
        for(i = 0; i < n; i++)
        {
	    delta.eye();
            s = S.row(i);
            for(j = 0; j < np; j++)
            {
                for(k = j + 1; k < np; k++)
                {
                    if( s(j) == s(k) )
                    {
                        delta(j,k) = delta(k,j) = 1;
                    } /* end if */

                } /* end loop k */
            } /* end loop j */

            phat += delta;
        } /* end loop i*/
        phat *= (1/double(n));

	// compute least squares score
        for(i = 0; i < n; i++)
        {
	    delta.eye();
            s = S.row(i);
            for(j = 0; j < np; j++)
            {
                for(k = j + 1; k < np; k++)
                {
                    if( s(j) == s(k) )
                    {
                        delta(j,k) = delta(k,j) = 1;
                    } /* end if */

                } /* end loop k */
            } /* end loop j */
            diffsq = pow((delta - phat),2);
            score(i) = sum(sum(diffsq));
        } /* end loop i*/
        // find index of minimum score and return
        ordscore = sort_index(score);
        minscoreloc = ordscore(0);

	// compute client cluster membership for chosen cluster
        ord = sort_index( trans( S.row(minscoreloc) ) );
        int M = Num(minscoreloc,0).n_rows;
	int startrow = 0; int endrow = 0; int nclustl;
	bigS.set_size(M,1);

        for(l = 0; l < M; l++)
        {
            nclustl = Num(minscoreloc,0)(l);
            endrow = startrow + nclustl - 1;
            bigS(l,0) = ord.rows(startrow,endrow) + 1; /* adding 1 shifts the subject index by 1 for R */
            startrow += nclustl;
        }

        END_RCPP
    } /* end function lsqcluster */

    // double s2 = std::inner_product(y.begin(),y.end(),y.begin(),double());

