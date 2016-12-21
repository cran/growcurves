#include "mmC.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


SEXP mmCplusDP(SEXP yvec, SEXP Xmat, SEXP Zmat, SEXP Wcase, SEXP Wperson, SEXP Omega,
        SEXP omegaplusvec, SEXP groupsvec, SEXP personsvec, SEXP niterInt,
        SEXP nburnInt, SEXP nthinInt, SEXP ustrengthd, SEXP shapealph, SEXP ratebeta)
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
    NumericMatrix Wcr(Wcase);
    NumericMatrix Wp(Wperson);
    NumericMatrix Or(Omega);
    NumericVector opr(omegaplusvec);
    IntegerVector gr(groupsvec);
    IntegerVector pr(personsvec);
    int niter = as<int>(niterInt);
    int nburn = as<int>(nburnInt);
    int nthin = as<int>(nthinInt);
    double ac = as<double>(shapealph);
    double bc = as<double>(ratebeta);
    double ustrength = as<double>(ustrengthd);
    int nkeep = (niter -  nburn)/nthin;


    // Extract key row and column dimensions we will need to sample parameters
    int nc =  Xr.nrow(), nf = Xr.ncol(), nr = Zr.ncol(), ns = Or.ncol();
    int ng = gr.length(); int npcbt = Wp.nrow();

    // Compute np, number of unique clients
    // algorithm assumes cases clustered by unique client
    IntegerVector dpr = diff(pr);
    Col<int32_t> diffpersons(dpr.begin(), nc - 1, false);
    int np = accu(diffpersons) + 1;
    /* int np = accumulate(dpr.begin,dpr.end,1); */

    // compute ng, number of session groups
    IntegerVector dgr = diff(gr);
    Col<int32_t> diffgroups(dgr.begin(),ns - 1, false);
    ng = accu(diffgroups) + 1;

    // Create armadillo structures we will use for our computations
    // We set copy_aux_mem = false from the constructor, meaning we re-use
    // the memory to in which our object is stored
    // Xr.begin() is a pointer, presumably of type double*
    mat xmat(Xr.begin(), nc, nf, false);
    mat zmat(Zr.begin(), nc, nr, false);
    mat wcase(Wcr.begin(),nc, ns, false);
    mat wpers(Wp.begin(),npcbt, ns, false); /* mmwt */
    mat omega(Or.begin(), ns, ns, false);
    colvec y(yr.begin(), nc, false);
    colvec omegaplus(opr.begin(), ns, false);
    Col<int32_t> persons(pr.begin(), nc, false);
    // set up data objects indexed by clustering unit
    field<mat> zsplit(np,1); field<mat> ztzsplit(np,1);
    field<mat> xsplit(np,1); field<mat> wsplit(np,1);
    field<colvec> ysplit(np,1); field<colvec> ytilsplit(np,1);
    CovUnitSplitMM(zmat, xmat, wcase, y, zsplit, ztzsplit, xsplit, wsplit, ysplit, persons);
    mat xtx; prodMatOne(xmat,xtx); 
    colvec wstws(ns); dotMMvecs(wcase, wstws); /* sequential, by effect, sampling of MM effects */

    // Set random number generator state
    RNGScope scope; /* Rcpp */
    //arma_rng::set_seed_random(); /* arma */

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
    mat MM(nkeep,npcbt);
    mat U(nkeep,ns);
    mat Taub(nkeep,nr);
    colvec Tauu(nkeep);
    colvec Taue(nkeep);
    icolvec numM(nkeep);
    imat S(nkeep,np);
    field<icolvec> Num(nkeep,1);
    field<ucolvec> bigS;
    field<mat> optPartition(3,1); /* Hold empirical probability of co-clustering matrix, index of L-sq clustering scores objects, and cluster identifiers per iteration */

    // Initialize parameter values
    /* cluster capture variables */
    double tauu, taubeta, taualph, taue;
    tauu = taubeta = taualph = taue = 1;
    rowvec taub(nr); taub.ones();
    IntegerVector ss = seq_len(np);
    icolvec s(ss.begin(), np, false);
    s -= 1; s = shuffle(s);
    int M = np; /* each person initialized with their own cluster */
    icolvec num(M); num.ones();
    mat    bstarmat = randn<mat>(M,nr)*sqrt(1/taub(0));
    mat    bmat(np,nr); bmat.zeros();
    colvec zb(nc); zb.zeros();
    mat qmat = -omega; qmat.diag() = omegaplus;/* prior precision for joint u */
    double conc = 1; /* DP concentration parameter */
    ucolvec ordscore; ordscore.zeros(); /* LS score order */
    mat phat(np,np); phat.zeros(); /* empirical co-clustering probability matrix */
    /* remaining parameter set */
    colvec beta = randn<colvec>(nf)*sqrt(1/taubeta);
    colvec u = randn<colvec>(ns)*sqrt(1/tauu);
    colvec resid(nc); resid.zeros(); /* regression residual */
    colvec mm(npcbt); mm.zeros(); /* client effect mapped from session eff */
    mat pbmat(nr,nr); pbmat.zeros(); /* prior precision matrix for B */
    double alpha = rnorm( 1, 0, sqrt(1/taualph) )[0];
    /* goodness of fit related measures */
    double deviance; colvec devres(4); rowvec devmarg(nc);
    rowvec logcpo(nc); logcpo.zeros(); double lpml;
    /* capture samples to return */
    int oo = 0, kk;
    
    // set hyperparameter values
    double a1, a2, a4, a7, b1, b2, b4, b7;
    a1 = b1 = ustrength; /* tauu */
    a2 = b2 = 1; /* taub */
    a4 = b4 = 1; /* taue */
    a7 = ac; b7 = bc; /* c */

    // Conduct posterior samples and store results
    int k, l;
    for(k = 0; k < niter; k++)
    {
        //if( (k % 1000) == 0 ) cout << "Interation: " << k << endl;
        clustermmstep(zsplit, ztzsplit, xsplit, wsplit, ysplit, ytilsplit, pbmat, 
                u, beta, alpha,
                taue, taub, bstarmat, s, num, M, conc,
                nr);
        concstep(conc, M, np, a7, b7);
        bstarstep(zsplit, ztzsplit, ytilsplit, taue, pbmat, s, num, bstarmat,
                zb, bmat, np, nr);
        betalscholstep(xmat, xtx, wcase, y, beta, u, alpha, taue, zb, nf);
        ustep(xmat, omega, wcase, wstws, wpers, beta, zb, y, omegaplus, u, mm,
                alpha, taue, tauu, ns, nc);
        alphalsstep(xmat, wcase, beta, zb, y, u, resid, alpha, taue, nc);
        taummdpstep(bstarmat, resid, qmat, u, tauu, taub, taue,
                a1, a2, a4, b1, b2, b4, M, ns, nc, ng);
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
                U.row(oo) = trans(u);
                MM.row(oo) = trans(mm);
                Tauu(oo) = tauu;
                Taue(oo) = taue;
                Taub.row(oo) = taub;
                Resid.row(oo) = trans(resid);
                numM(oo) = M;
                S.row(oo) = trans(s);
                Num(oo,0) = num;
                oo += 1; /* increment sample counter */
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
    cpo(Devmarg, logcpo, lpml);

    // Return results
    return Rcpp::List::create(Rcpp::Named("Deviance") = Deviance,
                                  Rcpp::Named("devres") = devres,
                                  //Rcpp::Named("logcpo") = logcpo,
                                  Rcpp::Named("lpml") = lpml,
				  Rcpp::Named("Beta") = Beta,
				  Rcpp::Named("Alpha") = Alpha,
                                  Rcpp::Named("C") = Conc,
                                  Rcpp::Named("B") = B,
                                  Rcpp::Named("U") = U,
                                  Rcpp::Named("MM") = MM,
                                  Rcpp::Named("Tauu") = Tauu,
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
    SEXP clustermmstep(const field<mat>& zsplit, const field<mat>& ztzsplit,
            const field<mat>& xsplit, const field<mat>& wsplit, const field<colvec>& ysplit, field<colvec>& ytilsplit, mat& pbmat, 
            const colvec& u,
            const colvec& beta, double alpha, double taue, const rowvec& taub,
            mat& bstarmat,
            icolvec& s, icolvec& num, int& M, double& conc,
            int nr)
    {
        BEGIN_RCPP
        // sample cluster assignments, s(1), ..., s(np)
        // fixing all other parameters, including bstarmat
        mat zj, zjtzj, wj, phib(nr,nr);
        colvec cj, cstarj, yj, ytildej, bstarj(nr), eb(nr), hb(nr);
	int np = zsplit.n_rows;
        double logq0, q0, sweights;
        int j, l, m;
        for(j = 0; j < np; j++)
        {
            /* subtracting cj eliminates need to t-fer xj to sample bstar */
	    zj = zsplit(j,0); yj = ysplit(j,0); zjtzj = ztzsplit(j,0);
	    wj = wsplit(j,0);
            cj = alpha + xsplit(j,0)*beta + wj*u;
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
                logq0 += dnorm(zro, 0.0, sqrt(1/taub(i)), true)[0];
            }
            eb = taue*trans(zj)*ytildej;
            phib = taue*zjtzj + pbmat; 
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
                rmvnchol(zj, zjtzj, pbmat, yj, cj, bstarj, nr, taue);
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

        } /* end loop j for cluster assignment */
        END_RCPP
    } /* end function bstep for cluster assignments, s, and computing zb */


    SEXP taummdpstep(const mat& bstarmat, const colvec& resid, const mat& qmat,
            const colvec& u, double& tauu, rowvec& taub, double& taue,
            double a1,double a2, double a4, double b1, double b2, double b4,  
            int M, int ns, int nc, int ng)
    {
        BEGIN_RCPP
        double a, b;
        /* tauu */
        a = 0.5*(double(ns) - double(ng)) + a1;
        b = 0.5*( as_scalar(trans(u)*qmat*u) ) + b1;
        tauu = rgamma(1, a, (1/b))[0];
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
        a = 0.5*double(nc) + a4;
        b = 0.5*dot( resid, resid ) + b4;
        taue = rgamma(1, a, (1/b))[0];
        END_RCPP
    } /* end function taustep to sample precision parameters */

    SEXP dic8comp(const colvec& Deviance, const colvec& Alpha,
            const mat& Beta, const mat& B, const colvec& Taue,
            const mat& U, const mat& xmat, const mat& zmat, const mat& wcase,
            const colvec& y, icolvec& persons, colvec& devres, int np)
    {
        BEGIN_RCPP
        // compute MAP parameter estimates
        double dbar = mean(Deviance);
        double alpha = mean(Alpha);
        double taue = mean(Taue);
        rowvec beta = mean(Beta,0);
        rowvec u    = mean(U,0);
        // compute zbhat
        int nc = xmat.n_rows;
        int nkeep = B.n_rows;
        colvec zb(nc); zb.zeros();
        colvec resid(nc); rowvec dhatrow;
        colvec residp = y - alpha - xmat*trans(beta) - wcase*trans(u);
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

SEXP CovUnitSplitMM(const mat& zmat, const mat& xmat, const mat& wmat, const colvec& y, field<mat>& zsplit, field<mat>& ztzsplit, 
			field<mat>& xsplit, field<mat>& wsplit, field<colvec>& ysplit, icolvec& persons)
   {
	BEGIN_RCPP
	// initialize objects
	int np = zsplit.n_rows; 
	int startrow, endrow, nj, j;
	mat zj, xj, wj; colvec yj;

	// compute number of repeated measures for each unit.
 	// assumes repeated units listed in nested / block fashion under unit
	icolvec positions(np); positions.zeros();
        icolvec::iterator aa = persons.begin();
        icolvec::iterator bb = persons.end();
        for(icolvec::iterator i = aa; i != bb ; i++)
        {
           positions(*i-1) += 1;
        }

	startrow 	= 0;
    	for(j = 0; j < np; j++)
        {
            // store by-person data cuts for later sampling of clust locs
            /* extract by-person matrix views from data*/
            nj = positions(j);
            endrow = startrow + nj - 1;
            zj = zmat.rows(startrow,endrow);
            xj = xmat.rows(startrow,endrow);
	    wj = wmat.rows(startrow,endrow);
            yj = y.rows(startrow,endrow);
            zsplit(j,0) = zj; /* np rows */
	    ztzsplit(j,0) = zj.t()*zj;
	    xsplit(j,0) = xj;
	    wsplit(j,0) = wj;
	    ysplit(j,0) = yj;
	    startrow	+= nj;
	}
	END_RCPP
    } /* end function splitting design objects by clustering unit */
