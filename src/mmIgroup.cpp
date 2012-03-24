#include "mmC.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


SEXP mmIgroup(SEXP yvec, SEXP Xmat, SEXP Zmat, SEXP Wcase, SEXP Wperson, SEXP Mmat,
        SEXP personsvec, SEXP niterInt, SEXP nburnInt, SEXP ustrengthd)
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
    NumericMatrix Mr(Mmat);
    NumericMatrix Wcr(Wcase);
    NumericMatrix Wp(Wperson);
    IntegerVector pr(personsvec);
    int niter = as<int>(niterInt);
    int nburn = as<int>(nburnInt);
    double ustrength = as<double>(ustrengthd);
    int nkeep = niter -  nburn;

    // Extract key row and column dimensions we will need to sample parameters
    int nc =  Xr.nrow(), nf = Xr.ncol(), nr = Zr.ncol(), ns = Wcr.ncol();
    int npcbt = Wp.nrow(), ng = Mr.ncol();

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
    mat mmat(Mr.begin(), ns, ng, false);
    mat wcase(Wcr.begin(),nc, ns, false);
    mat wpers(Wp.begin(),npcbt, ns, false); /* mmwt */
    colvec y(yr.begin(), nc, false);
    icolvec persons(pr.begin(), nc, false);

    // Set random number generator state
    RNGScope scope; /* Rcpp */
    srand ( time(NULL) ); /* arma */

    // Armadillo structures to capture results
    // Will be wrapped into a list element on return
    colvec Deviance(nkeep); Deviance.zeros();
    mat Devmarg(nkeep,nc);
    mat Resid(nkeep,nc);
    colvec Alpha(nkeep);
    mat Beta(nkeep,nf);
    mat B(nkeep,nr*np);
    rowvec brow(np*nr); brow.zeros();
    mat MM(nkeep,npcbt);
    mat Eta(nkeep,ng);
    mat U(nkeep,ns);
    mat Taub(nkeep,nr);
    colvec Tauu(nkeep);
    colvec Taubeta(nkeep);
    colvec Taualph(nkeep);
    colvec Taue(nkeep);
    colvec Taueta(nkeep);

    // Initialize parameter values
    double tauu, taubeta, taualph, taue, taueta;
    tauu = taubeta = taualph = taue = taueta = 1;
    rowvec taub(nr); taub.ones();
    colvec beta = randn<colvec>(nf)*sqrt(1/taubeta);
    mat    bmat = randn<mat>(np,nr)*sqrt(1/taub(0));
    colvec eta(ng); eta.zeros(); /* group means of session effects */
    colvec u = randn<colvec>(ns)*sqrt(1/tauu);
    colvec zb(nc); zb.zeros();
    colvec resid(nc); resid.zeros(); /* regression residual */
    colvec mm(npcbt); mm.zeros(); /* client effect mapped from session eff */
    double alpha = rnorm( 1, 0, sqrt(1/taualph) )[0];
    /* goodness of fit related measures */
    double deviance; colvec devres(4); rowvec devmarg(nc);
    rowvec logcpo(nc); logcpo.zeros(); double lpml;

    // set hyperparameter values
    double a1, a2, a4, a5, a6, a7, b1, b2, b4, b5, b6, b7;
    a1 = b1 = ustrength; /* tauu */
    a7 = b7 = 1; /* taueta*/
    a2 = b2 = 1; /* taub */
    a4 = b4 = 1; /* taue */
    a5 = 0.01; b5 = 0.01; /* taualph */
    a6 = 0.1; b6 = 0.1; /* taubeta */

    // Conduct posterior samples and store results
    int k, l;
    for(k = 0; k < niter; k++)
    {
        //if( (k % 1000) == 0 ) cout << "Interation: " << k << endl;
        bcholstep(xmat, zmat, wcase, y, beta, u, alpha, taue, taub,
            persons, zb, bmat, np, nr);
        betacholstep(xmat, wcase, y, beta, u, alpha, taue, taubeta, zb, nf);
        uindetastep(xmat, wcase, wpers, mmat, beta, zb, y, eta, u, mm,
                alpha, taue, tauu, ns);
        etastep(mmat, u, eta, tauu, taueta);
        alphastep(xmat, wcase, beta, zb, y, u, resid, alpha, taue, taualph, nc);
        tauetastep(bmat, mmat, resid, u, eta, beta, tauu, taub, taue,
                taualph, taubeta, taueta, alpha, a1, a2, a4, a5, a6, a7, b1, b2,
                b4, b5, b6, b7, ns, np, nc, nf, ng);
        if(k >= nburn)
        {
            deviance =  dev(resid, taue);
            dmarg(resid, taue, devmarg); /* 1 x nc vector of densities*/
            int oo = k - nburn;
            Deviance(oo) = deviance;
            Devmarg.row(oo)  = devmarg;
            Beta.row(oo) = trans(beta);
            Alpha(oo) = alpha;
            for( l = 0; l < nr; l++)
            {
                brow.cols( np*l,(np*(l+1)-1) ) = trans( bmat.col(l) );
            }
            B.row(oo) = brow;
            U.row(oo) = trans(u);
            Eta.row(oo) = trans(eta);
            MM.row(oo) = trans(mm);
            Tauu(oo) = tauu;
            Taubeta(oo) = taubeta;
            Taualph(oo) = taualph;
            Taue(oo) = taue;
            Taueta(oo) = taueta;
            Taub.row(oo) = taub;
            Resid.row(oo) = trans(resid);
        } /* end if k > burnin, record results for return */

    } /* end MCMC loop over k */

    // DIC
    dic1comp(Deviance, Alpha, Beta, B, Taue, U, xmat, zmat, wcase,
            y, persons, devres, np);
    cpo(Devmarg, logcpo, lpml);

    // Return results
    return Rcpp::List::create(Rcpp::Named("Deviance") = Deviance,
                                  Rcpp::Named("devres") = devres,
                                  Rcpp::Named("logcpo") = logcpo,
                                  Rcpp::Named("lpml") = lpml,
				  Rcpp::Named("Beta") = Beta,
				  Rcpp::Named("Alpha") = Alpha,
                                  Rcpp::Named("B") = B,
                                  Rcpp::Named("U") = U,
                                  Rcpp::Named("Eta") = Eta,
                                  Rcpp::Named("MM") = MM,
                                  Rcpp::Named("Tauu") = Tauu,
                                  Rcpp::Named("Taubeta") = Taubeta,
                                  Rcpp::Named("Taualph") = Taualph,
                                  Rcpp::Named("Taue") = Taue,
                                  Rcpp::Named("Taueta") = Taueta,
				  Rcpp::Named("Taub") = Taub,
                                  Rcpp::Named("Residuals") = Resid
				  );

    // Print CPU runtime
     //double n_secs = timer.toc();
     //cout << "took " << n_secs << " seconds" << endl;
END_RCPP
} /* end MCMC function returning SEXP */

SEXP uindetastep(const mat& xmat, const mat& wcase, const mat& wpers, 
            const mat& mmat, const colvec& beta, const colvec& zb,
            const colvec& y, const colvec& eta, colvec& u,
            colvec& mm,  double alpha, double taue,
            double tauu, int ns)
    {
        BEGIN_RCPP
        // Sample posterior of u (ns x 1) from posterior Gaussian
        // compute prior precision matrix
        mat pu(ns,ns); pu.eye(); pu *= tauu;
        colvec c = alpha + xmat*beta + zb; /* nc x 1 */
        colvec m = mmat*eta; /* ns x 1 */
        rmvnmean(wcase, pu, y, c, m, u, ns, taue);

       // compute npcbt x 1, mm, resulting from mapping u to the npcbt clients
        mm = wpers*u;

        END_RCPP
    } /* end function to sample nf x 1 fixed effects, beta */

SEXP etastep(const mat& mmat, 
            const colvec& u, colvec &eta,
            double tauu, double taueta)
    {
        BEGIN_RCPP
        // Sample posterior of u (ns x 1) from posterior Gaussian
        // compute prior precision matrix
        int ns = u.n_elem, ng = eta.n_elem;
        mat pu(ng,ng); pu.eye(); pu *= taueta;
        colvec c(ns); c.zeros();
        rmvnchol(mmat, pu, u, c, eta, ng, tauu);
        // enforce constraint for mean(eta) = 0
        eta -= mean(eta);

        END_RCPP
    } /* end function to sample nf x 1 fixed effects, beta */

SEXP tauetastep(const mat& bmat, const mat& mmat,
            const colvec& resid, const colvec& u,
            const colvec& eta, const colvec& beta, double& tauu, rowvec& taub,
            double& taue, double& taualph, double& taubeta, double& taueta,
            double alpha,
            double a1, double a2, double a4, double a5, double a6, double a7,
            double b1, double b2,  double b4, double b5, double b6, double b7,
            int ns, int np, int nc, int nf, int ng)
    {
        BEGIN_RCPP
        double a, b;
        /* tauu */
        a = 0.5*double(ns) + a1;
        colvec uresid = u - mmat*eta;
        b = 0.5* dot( uresid , uresid ) + b1;
        tauu = rgamma(1, a, (1/b))[0];
        /* taueta */
        a = 0.5*double(ng) + a7;
        b = 0.5* dot( eta , eta ) + b7;
        taueta = rgamma(1, a, (1/b))[0];
        /* tau0 */
        int i;
        int nr = bmat.n_cols;
        a = 0.5*double(np) + a2;
        for(i = 0; i < nr; i++)
        {
            b =  0.5*dot( bmat.col(i), bmat.col(i) ) + b2;
            taub(i) = rgamma(1, a, (1/b))[0];
        }
        /* tau1 */
        //a = 0.5*double(np) + a3;
        // mat bless0mat = bmat.cols(1,nr-1);
        //b = 0.5*as_scalar( sum( diagvec(bless0mat*trans(bless0mat)) ) ) + b3;
        //tau1 = rgamma(1, a, (1/b))[0];
        /* taue */
        a = 0.5*double(nc) + a4;
        b = 0.5*dot( resid, resid ) + b4;
        taue = rgamma(1, a, (1/b))[0];
        /* taualph */
        a = 0.5 + a5;
        b = 0.5*pow(alpha,2) + b5;
        taualph = rgamma(1, a, (1/b))[0];
        /* taubeta */
        a = 0.5*double(nf) + a6;
        b = 0.5*dot(beta,beta) + b6;
        taubeta = rgamma(1, a, (1/b))[0];
        END_RCPP
    } /* end function taustep to sample precision parameters */
