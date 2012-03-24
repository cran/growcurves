#include "mmC.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


SEXP lgm(SEXP yvec, SEXP Xmat, SEXP Zmat, SEXP personsvec, SEXP niterInt,
         SEXP nburnInt, SEXP nthinInt)
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

    // Set random number generator state
    RNGScope scope; /* Rcpp */
    srand ( time(NULL) ); /* arma */

    // Armadillo structures to capture results
    // Will be wrapped into a list element on return
    colvec Deviance(nkeep); Deviance.zeros();
    colvec Alpha(nkeep);
    mat Devmarg(nkeep,nc);
    mat Resid(nkeep,nc);
    mat Beta(nkeep,nf);
    mat B(nkeep,nr*np);
    rowvec brow(np*nr); brow.zeros();
    mat Taub(nkeep,nr);
    colvec Taue(nkeep);

    // Initialize parameter values
    double taubeta, taualph, taue;
    taubeta = taualph = taue = 1;
    rowvec taub(nr); taub.ones();
    colvec beta = randn<colvec>(nf)*sqrt(1/taubeta);
    mat    bmat = randn<mat>(np,nr)*sqrt(1/taub(0));
    colvec zb(nc); zb.zeros();
    colvec resid(nc); resid.zeros(); /* regression residual */
    double alpha = rnorm( 1, 0, sqrt(1/taualph) )[0];
    /* goodness of fit related measures */
    double deviance = 0; colvec devres(4); devres.zeros();
    colvec devres3(4); devres3.zeros();
    rowvec devmarg(nc);
    rowvec logcpo(nc); logcpo.zeros(); double lpml;
    /* capture samples to return */
    int oo = 0, kk;

    // set hyperparameter values
    double a2, a4, b2, b4;
    a2 = b2 = 1; /* taub */
    a4 = b4 = 1; /* taue */

    // Conduct posterior samples and store results
    int k, l;
    for(k = 0; k < niter; k++)
    {
        //if( (k % 1000) == 0 ) cout << "Interation: " << k << endl;
        blgmstep(xmat, zmat, y, beta, alpha, taue, taub,
            persons, zb, bmat, np, nr);
        betalgmstep(xmat, y, beta, alpha, taue, zb, nf);
        alphalgmstep(xmat, beta, zb, y, resid, alpha, taue, nc);
        taulgmstep(bmat, resid, taub, taue, a2, a4, b2, b4, np, nc);
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
                for( l = 0; l < nr; l++)
                {
                    brow.cols( np*l,(np*(l+1)-1) ) = trans( bmat.col(l) );
                }
                B.row(oo) = brow;
                Taue(oo) = taue;
                Taub.row(oo) = taub;
                Resid.row(oo) = trans(resid);
                oo += 1; /* increment sample return counter */
            } /* end conditional statement on returning sample */
        } /* end if k > burnin, record results for return */

    } /* end MCMC loop over k */

    // DIC
    dic1lgmcomp(Deviance, Alpha, Beta, B, Taue, xmat, zmat,
            y, persons, devres, np); /* devres = c(dic,dbar,dhat,pd) */
    dic3comp(Deviance, Devmarg, devres3);
    cpo(Devmarg, logcpo, lpml);

    // Return results
    return Rcpp::List::create(Rcpp::Named("Deviance") = Deviance,
                                  Rcpp::Named("devres") = devres,
                                  Rcpp::Named("devres3") = devres3,
                                  Rcpp::Named("logcpo") = logcpo,
                                  Rcpp::Named("lpml") = lpml,
				  Rcpp::Named("Beta") = Beta,
				  Rcpp::Named("Alpha") = Alpha,
                                  Rcpp::Named("B") = B,
                                  Rcpp::Named("Taue") = Taue,
                                  Rcpp::Named("Taub") = Taub,
                                  Rcpp::Named("Residuals") = Resid
				  );
    // Print CPU runtime
     //double n_secs = timer.toc();
     //cout << "took " << n_secs << " seconds" << endl;
END_RCPP
} /* end MCMC function returning SEXP */

 SEXP blgmstep(const mat& xmat, const mat& zmat,
            const colvec& y, const colvec& beta,
            double alpha, double taue,
            const rowvec& taub, icolvec& persons, colvec& zb,
            mat& bmat, int np, int nr)
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

        // Compute cholesky of prior precision, pbmat, for mvn sampling
        // of bmat through cholesky decomposition
        mat pbmat; pbmat.eye(nr,nr);
        int i;
        for(i = 0; i < nr; i++)
        {
            pbmat.diag()(i) = taub(i);
        }
        mat zj, xj, yj;
        colvec cj, bj(nr);

        // sample bmat, independently, by person
        int startrow = 0; int persrow = 0;
        int j, nj, endrow;
        for(j=0; j < np; j++)
        {
            /* extract by-person matrix views from data*/
            nj = positions(j);
            endrow = startrow + nj - 1;
            zj = zmat.rows(startrow,endrow);
            xj = xmat.rows(startrow,endrow);
            yj = y.rows(startrow,endrow);

            /* sample posterior of nj block of bmat*/
            cj = alpha + xj*beta;
            bj.zeros();
            rmvnchol(zj, pbmat, yj, cj, bj, nr, taue);
            bmat.row(persrow) = trans(bj);

            // compute zbj (jth piece) of zb = [z1*b1vec,...,znp*bnpvec]
            zb.rows(startrow,endrow) = zj*bj;

            //re-set start positions
            startrow += nj;
            persrow++; /*bmat has np rows, increment once for each iteration*/
        } /* end loop j for sampling bj*/
        END_RCPP
    } /* end function bstep for sampling bmat and zb */

    SEXP betalgmstep(const mat& xmat, const colvec& y,
                colvec& beta, double alpha, double taue,
                const colvec& zb, int nf)
    {
        BEGIN_RCPP
        // Sample posterior of beta (nf x 1) from posterior Gaussian
        colvec c = alpha + zb; /* nc x 1 */
        rmvnlschol(xmat, y, c, beta, nf, taue);
        END_RCPP
    } /* end function to sample nf x 1 fixed effects, beta */

    SEXP alphalgmstep(const mat& xmat, const colvec& beta,
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

    SEXP taulgmstep(const mat& bmat, const colvec& resid, 
            rowvec& taub, double& taue, double a2, double a4, 
            double b2, double b4, int np, int nc)
    {
        BEGIN_RCPP
        double a, b;
        /* taub */
        int nr = bmat.n_cols;
        a = 0.5*double(np) + a2;
        int i;
        for(i = 0; i < nr; i++)
        {
            b =  0.5*dot( bmat.col(i), bmat.col(i) ) + b2;
            taub(i) = rgamma(1, a, (1/b))[0];
        }
        /* taue */
        a = 0.5*double(nc) + a4;
        b = 0.5*dot( resid, resid ) + b4;
        taue = rgamma(1, a, (1/b))[0];
        END_RCPP
    } /* end function taustep to sample precision parameters */


    SEXP dic1lgmcomp(const colvec& Deviance, const colvec& Alpha,
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
        rowvec b = mean(B,0);
        // compute zbhat
        int nc = xmat.n_rows;
        colvec zb(nc); zb.zeros();
        zbcomp(b, zmat, persons, zb, np);
        colvec resid = y - alpha - xmat*trans(beta) - zb;
        double dhat = dev(resid,taue);
        double dic = 2*dbar - dhat;
        double pd = dic - dbar;
        devres(0) = dic; devres(1) = dbar; devres(2) = dhat; devres(3) = pd;
        END_RCPP
    }
