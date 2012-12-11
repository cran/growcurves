#include "mmC.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


SEXP mmCmvplusDP(SEXP yvec, SEXP Xmat, SEXP Zmat, SEXP Hmat, SEXP Wcase, SEXP Wperson,
        SEXP Omega, SEXP omegaplusvec, SEXP groupsvec, SEXP personsvec,
        SEXP niterInt, SEXP nburnInt, SEXP nthinInt, SEXP ustrengthd, SEXP corsessInt,
        SEXP shapealph, SEXP ratebeta, SEXP typeMM)
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
    NumericMatrix Hr(Hmat);
    NumericMatrix Wcr(Wcase);
    NumericMatrix Wp(Wperson);
    NumericMatrix Or(Omega);
    NumericVector opr(omegaplusvec);
    IntegerVector gr(groupsvec);
    IntegerVector pr(personsvec);
    int niter = as<int>(niterInt);
    int nburn = as<int>(nburnInt);
    int nthin = as<int>(nthinInt);
    int corsess = as<int>(corsessInt); /* correlation between sets of sess eff*/
    int typemm = as<int>(typeMM); /* typemm == 0 for ind prior on umat, 1 for CAR prior on umat */
    double ac = as<double>(shapealph);
    double bc = as<double>(ratebeta);
    double ustrength = as<double>(ustrengthd);
    int nkeep = (niter -  nburn)/nthin;

    // Extract key row and column dimensions we will need to sample parameters
    int nc = Xr.nrow(), nf = Xr.ncol(), nr = Zr.ncol(), nv = Hr.ncol(), ns = Wcr.ncol();
    int npcbt = Wp.nrow(); int nrho = 0.5*nv*(nv-1);

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
    mat hmat(Hr.begin(), nc, nv, false);
    mat wcase(Wcr.begin(),nc, ns, false);
    mat wpers(Wp.begin(),npcbt, ns, false); /* mmwt */
    mat omega(Or.begin(), ns, ns, false);
    colvec y(yr.begin(), nc, false);
    colvec omegaplus(opr.begin(), ns, false);
    // icolvec groups(gr.begin(), ns, false);
    icolvec persons(pr.begin(), nc, false);
    // set up data objects indexed by clustering unit
    field<mat> zsplit(np,1); field<mat> ztzsplit(np,1);
    field<mat> xsplit(np,1); field<mat> wsplit(np,1);
    field<mat> hsplit(np,1);
    field<colvec> ysplit(np,1); field<colvec> ytilsplit(np,1);
    CovUnitSplitMV(zmat, xmat, wcase, hmat, y, zsplit, ztzsplit, xsplit, wsplit, hsplit, ysplit, persons);
    mat xtx; prodMatOne(xmat,xtx); 
    field<mat> hmatws(ns,1), hwsthws(ns,1);
    prodHW(wcase, hmat, hmatws, hwsthws); /* multiple each nc x nv, hmat matrix by column s of wcase, ws - data object to use for sampling MM effects  */

    // Set random number generator state
    RNGScope scope; /* Rcpp */
    srand ( time(NULL) ); /* arma */

    // Armadillo structures to capture results
    // Will be wrapped into a list element on return
    colvec Deviance(nkeep); Deviance.zeros();
    mat Devmarg(nkeep,nc);
    mat Resid(nkeep,nc);
    colvec Alpha(nkeep);
    colvec Conc(nkeep);
    mat Beta(nkeep,nf);
    mat B(nkeep,nr*np);
    rowvec brow(np*nr); brow.zeros();
    mat MM(nkeep,nv*npcbt);
    rowvec mmrow(npcbt*nv); mmrow.zeros();
    mat U(nkeep,nv*ns);
    rowvec urow(ns*nv); urow.zeros();
    mat Taub(nkeep,nr);
    mat Tauu(nkeep,nv);
    rowvec taurow(nv); taurow.zeros();
    mat Rhotauu(nkeep,nrho);
    rowvec rhotaurow(nrho); rhotaurow.zeros();
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
    ucolvec ordscore; ordscore.zeros(); /* LS score order */
    mat phat(np,np); phat.zeros(); /* empirical co-clustering probability matrix */
    // capture session effects in ns x nsr matrix, umat
    mat Sigma(nv,nv); /* variance of u */
    mat L(nv,nv); int i, j; /* precision for umat */
    for(i = 0; i < nv; i++)
    {
        Sigma(i,i) = 1/ustrength;
    }
    for(i = 0; i < nv; i++)
    {
        for( j = i+1; j < nv; j++)
        {
            Sigma(i,j) = Sigma(j,i) = corsess*sqrt(Sigma(i,i))*sqrt(Sigma(j,j));
        } /*end loop j */
    } /* end loop i filling variance matrix for umat, S */
    L = inv(Sigma); mat V(nv,nv); V = L; /* V is prior for L */
    double nu = nv + 1; /* df for L ~ W(V,nu) */
    /* compute correlation matrix for each iteration */
    mat Linv(nv,nv); int nelem, totelem; double rho;
    mat umat(ns,nv); rowvec m(nv); m.zeros(); rmvnrnd(L,umat,m);
    mat mmmat(npcbt,nv); mmmat.zeros(); /* client effect mapped from session */
    int ng; /* number of groups for MM (session) effects */
    mat qmat;
    if( typemm == 0)  /* mmi */
    {
	qmat.eye(ns,ns);
	ng = 0;
    }else{ /* typemm = 1 */
    	qmat = -omega; qmat.diag() = omegaplus;/* prior precision for joint u */
	// compute ng, number of session groups
    	IntegerVector dgr = diff(gr);
    	icolvec diffgroups(dgr.begin(),ns - 1, false);
    	ng = accu(diffgroups) + 1;
    }
    /* remaining parameters */
    colvec beta = randn<colvec>(nf)*sqrt(1/taubeta);
    colvec resid(nc); resid.zeros(); /* regression residual */
    double alpha = rnorm( 1, 0, sqrt(1/taualph) )[0];
    mat pbmat(nr,nr); pbmat.zeros();
    /* goodness of fit related measures */
    double deviance; colvec devres(4); rowvec devmarg(nc);
    rowvec logcpo(nc); logcpo.zeros(); double lpml;
    /* capture samples to return */
    int oo = 0, kk;

    // set hyperparameter values
    double a2, a4, a7, b2, b4, b7;
    a2 = b2 = 1; /* taub */
    a4 = b4 = 1; /* taue */
    a7 = ac; b7 = bc; /* conc */

    // Conduct posterior samples and store results
    int k,l,h;
    for(k = 0; k < niter; k++)
    {
        //if( (k % 1000) == 0 ) Rcout << "Interation: " << k << endl;
        clustermvstep(zsplit, ztzsplit, xsplit, wsplit, hsplit, ysplit, 
		ytilsplit, pbmat, umat, beta, alpha, taue, taub, 
		bstarmat, s, num, M, conc, nr);
        concstep(conc, M, np, a7, b7);
        bstarstep(zsplit, ztzsplit, ytilsplit, taue, pbmat, s, num, bstarmat,
                zb, bmat, np, nr);
        betamvlsstep(xmat, xtx, wcase, hmat, y, beta, umat, alpha, taue, 
                    zb, nf);
	if( typemm == 0)
	{
		uindmvstep(xmat, wcase, wpers, hmat, hmatws, hwsthws, zb, beta, y, umat,
                	mmmat, alpha, taue, L, ns, nc);
	}else{ /* typemm == 1, MM-MCAR */
		umvstep(xmat, omega, wcase, wpers, hmat, hmatws, hwsthws, zb, beta, y, omegaplus, umat,
                	mmmat, alpha, taue, L, ns, nc);
	}
        alphamvlsstep(xmat, wcase, beta, hmat, zb, y, umat, resid, alpha,
                    taue, nc);
        taumvdpstep(bstarmat, resid, qmat, umat, V, L, taub, taue,
                nu, a2, a4, b2, b4, M, ns, nc, ng);
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
            	/* nv sets of session effects */
            	for( h = 0; h < nv; h++)
            	{
                	urow.cols( ns*h,(ns*(h+1)-1) ) = trans( umat.col(h) );
                	mmrow.cols( npcbt*h,(npcbt*(h+1)-1) ) = trans( mmmat.col(h) );
                	taurow.col(h) = L(h,h);
            	}
            	/* correlation coefficients between sets of session effects */
            	rho = 0; Linv = inv(L); totelem = 0; nelem = 0;
            	for( h = 0; h < (nv-1); h++)
            	{
                	nelem = nv - (h+1);
                	for(l = (h+1); l < nv; l++)
                	{
                    		rho = Linv(h,l)/( sqrt(Linv(h,h))*sqrt(Linv(l,l)) );
                    		rhotaurow.col(totelem + (l - h) - 1) = rho;
                	}
                	totelem += nelem;
            	}
            	B.row(oo) = brow;
            	U.row(oo) = urow;
            	MM.row(oo) = mmrow;
            	Tauu.row(oo) = taurow;
            	Rhotauu.row(oo) = rhotaurow;
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
                                  // Rcpp::Named("logcpo") = logcpo,
                                  Rcpp::Named("lpml") = lpml,
				  Rcpp::Named("Beta") = Beta,
				  Rcpp::Named("Alpha") = Alpha,
                                  Rcpp::Named("C") = Conc,
                                  Rcpp::Named("B") = B,
                                  Rcpp::Named("U") = U,
                                  Rcpp::Named("MM") = MM,
                                  Rcpp::Named("Tauu") = Tauu,
                                  Rcpp::Named("Rhotauu") = Rhotauu,
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
    SEXP clustermvstep(const field<mat>& zsplit, const field<mat>& ztzsplit,
	    const field<mat>& xsplit, const field<mat>& wsplit, const field<mat>& hsplit,
            const field<colvec>& ysplit, field<colvec>& ytilsplit, 
            mat& pbmat, const mat& umat,
            const colvec& beta, double alpha, double taue, const rowvec& taub,
            mat& bstarmat,
            icolvec& s, icolvec& num, int& M, double& conc,
            int nr)
    {
        BEGIN_RCPP
        // sample cluster assignments, s(1), ..., s(np)
        // fixing all other parameters, including bstarmat
        mat zj, zjtzj, hj, xj, wj, phib(nr,nr);
        colvec cj, cstarj, yj, ytildej, bstarj(nr), eb(nr), hb(nr);
	int np = zsplit.n_rows;
        double logq0, q0, sweights;
        int i, j, l, m;
        for(j = 0; j < np; j++)
        {
            // store by-person data cuts for later sampling of clust locs
            /* extract by-person matrix views from data*/
            zj = zsplit(j,0); yj = ysplit(j,0); zjtzj = ztzsplit(j,0);
	    wj = wsplit(j,0); hj = hsplit(j,0); xj = xsplit(j,0);

            /* subtracting cj eliminates need to t-fer xj to sample bstar */
            cj = alpha + xj*beta;
            int nv = umat.n_cols;
            for(i = 0; i < nv; i++)
            {
                cj += hj.col(i) % (wj*umat.col(i));
            }
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
            eb 		= taue*trans(zj)*ytildej;
            phib 	= taue*zjtzj + pbmat; 
	    hb 		= solve(phib,eb);
            logq0 	-= logdens(hb,phib); /* dmvn(hb,phib^-1) */
            q0 		= exp(logq0);

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

    SEXP taumvdpstep(const mat& bstarmat, const colvec& resid, const mat& qmat,
            const mat& umat, const mat& V, mat& L, rowvec& taub, double& taue,
            double nu, double a2, double a4, double b2, double b4, int M, int ns,
            int nc, int ng)
    {
        BEGIN_RCPP
        double a, b;
        /* L */
        mat Vstar = inv( trans(umat)*qmat*umat + inv(V) );
        double nustar = double(ns) - double(ng) + nu;
        wishrnd(L, Vstar, nustar);
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

    SEXP CovUnitSplitMV(const mat& zmat, const mat& xmat, const mat& wmat, const mat& hmat, const colvec& y, field<mat>& zsplit, field<mat>& ztzsplit, 
			field<mat>& xsplit, field<mat>& wsplit, field<mat>& hsplit, field<colvec>& ysplit, icolvec& persons)
   {
	BEGIN_RCPP
	// initialize objects
	int np = zsplit.n_rows; 
	int startrow, endrow, nj, j;
	mat zj, xj, wj, hj; colvec yj;

	// compute number of repeated measures for each unit.
 	// assumes repeated units listed in nested / block fashion under unit
	icolvec positions(np); positions.zeros();
        icolvec::iterator aa = persons.begin();
        icolvec::iterator bb = persons.end();
        for(icolvec::iterator i = aa; i != bb ; i++)
        {
           positions(*i-1) += 1;
        }

	startrow = 0;
    	for(j = 0; j < np; j++)
        {
            // store by-person data cuts for later sampling of clust locs
            /* extract by-person matrix views from data*/
            nj = positions(j);
            endrow = startrow + nj - 1;
            zj = zmat.rows(startrow,endrow);
            xj = xmat.rows(startrow,endrow);
	    wj = wmat.rows(startrow,endrow);
	    hj = hmat.rows(startrow,endrow);
            yj = y.rows(startrow,endrow);
            zsplit(j,0) = zj; /* np rows */
	    ztzsplit(j,0) = zj.t()*zj;
	    xsplit(j,0) = xj;
	    wsplit(j,0) = wj;
	    hsplit(j,0) = hj;
	    ysplit(j,0) = yj;
	    startrow	+= nj;
	}
	END_RCPP
    } /* end function splitting design objects by clustering unit */
