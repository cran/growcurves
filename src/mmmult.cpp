#include "mmC.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


SEXP mmmult(SEXP yvec, SEXP Xmat, SEXP Zmat, SEXP WcaseList, SEXP MmatList, SEXP OmegaList,
        SEXP omegapluslist, SEXP ngsvec, SEXP personsvec,  SEXP typeTreat, SEXP niterInt, SEXP nburnInt, 
        SEXP nthinInt, SEXP shapealph, SEXP ratebeta, SEXP ustrengthd)
{
BEGIN_RCPP
    NumericVector yr(yvec);
    IntegerVector numgs(ngsvec);
    NumericMatrix Xr(Xmat);
    NumericMatrix Zr(Zmat);
    List Wcase(WcaseList);
    List Mmat(MmatList);
    List Omega(OmegaList); /* List of Omega matrices for CAR priors */
    List omegaplus(omegapluslist);
    // vector of covariance form choices for each treatment type under a multivariate Gaussian.
    // 1 = CAR, 2 = MVN, 3 = IND
    IntegerVector typeT(typeTreat); /* same length as numT */
    IntegerVector pr(personsvec);
    int niter = as<int>(niterInt);
    int nburn = as<int>(nburnInt);
    int nthin = as<int>(nthinInt);
    int nkeep = (niter -  nburn)/nthin;
    int nOmega = Omega.size(); /* should equal numcar */
    double ac = as<double>(shapealph);
    double bc = as<double>(ratebeta);
    double ustrength = as<double>(ustrengthd);

    // Extract key row and column dimensions we will need to sample parameters
    int nc =  Xr.nrow(), nf = Xr.ncol(), nr = Zr.ncol(), nty = typeT.size();
    /* nty = number of treatment types */
    
    // Compute np, number of unique persons
    /* algorithm assumes cases clustered by unique client */ 
    IntegerVector dpr = diff(pr);
    icolvec diffpersons(dpr.begin(), nc - 1, false);
    int np = accu(diffpersons) + 1;
    
    // Create armadillo structures we will use for our computations
    // We set copy_aux_mem = false from the constructor, meaning we re-use
    // the memory to in which our object is stored
    // Xr.begin() is a pointer, presumably of type double*
    mat xmat(Xr.begin(), nc, nf, false);
    mat zmat(Zr.begin(), nc, nr, false);
    colvec y(yr.begin(), nc, false);
    icolvec persons(pr.begin(), nc, false);
    /* use for transferring by-subject objects from clustering to sampling locations for DP(subjects) */
    field<mat> zsplit(np,1);
    field<colvec> ytilsplit(np,1);
    /* objects for sampling sessions, us */
    icolvec typet(typeT.begin(), nty, false);
    field<mat> wmats(nty,1);  /* multiple membership matrices mapping sessions to subjects for each treatment */
    icolvec numt(nty);
    for(int i = 0; i < nty; i++)
    {
             wmats(i,0) = as<mat>(Wcase[i]);
             numt(i) = wmats(i,0).n_cols;
    } 
    /* ensure all the wmats have equal number of rows = nc */
    for(int i = 0; i < nty; i++)
    {
        for(int j = (i+1); j < nty; j++)
        {
            if(wmats(i,0).n_rows != wmats(j,0).n_rows)
            {
                throw std::range_error("n_rows of all wmats must be equal to number of cases");
            }
        }
    }
    
    // capture number of MM terms (treatments).
    // use to dimension return objects.
    uvec poscar = find(typet == 1); int numcar = poscar.n_elem;
    uvec posigp = find(typet == 3); int numigp = posigp.n_elem;
    uvec posdp  = find(typet == 4); int numdp  = posdp.n_elem;

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
    // session effect sampling and return structures 
    int i, dim;
    mat Tauu(nkeep,nty); /* return object for tauus for each typet */
    field<mat> U(nty,1); /* List element is typet(m).  nkeep x dim_m matrix of sampled session effects (dosages) */
    for(i = 0; i < nty; i++)
    {
        dim = numt(i);
        U(i,0).set_size(nkeep,dim);
    }
     
    // Session parameters

    /* sampling objects */
    field<colvec> us(nty,1); /* groups all MM objects */
    rowvec tauus(nty); tauus.ones(); /* groups all MM objects */
    field<colvec> etas; /* mmigrp */
    rowvec tauetas; /* mmigrp */
    field<mat> mmats; /* S x G session mapping to group for mmigrp */
    icolvec ngs;  /*  mmcar */
    field<mat> omegamats;  /* CAR */
    field<rowvec> omegapluses; /* CAR */
    field<colvec> ustars; /* dp */
    field<ucolvec> ss; /* dp */
    field<icolvec> nums; /* dp */
    irowvec Ms; /* dp */
    IntegerVector slabs; /* dp */
    ucolvec su; /* dp */
    rowvec concs; /* dp */
    
    /* return objects */
    imat numMu; /* dp */
    mat Concu; /* dp */
    field<umat> Sulist; /* dp */
    field<umat> bigSulist; /* dp */ 
    field<mat> Etalist; /* igrp */

    // set initial values for sampling objects and dimensions for return objects 
    /* session effects, us, initial values */
    for(i = 0; i < nty; i++)
    {
        int nsessi     = numt(i);
        us(i,0) = randn<colvec>(nsessi)*sqrt(10);
    }
    if(numdp > 0) /* return precision matrix objects for mvn structures */
    {
       //Define field object of dp parameters for sampling */
       /* sampling objects */
       concs.set_size(numdp); concs.ones();
       ustars.set_size(numdp,1); 
       ss.set_size(numdp,1);   
       nums.set_size(numdp,1);
       Ms.set_size(numdp);
       /* return objects */
       numMu.set_size(nkeep, numdp);
       Concu.set_size(nkeep, numdp);
       Sulist.set_size(numdp,1);
       bigSulist.set_size(numdp,1);
       // fill in each field member with a vector of tau values for that type */   
        for(i = 0; i < numdp; i++)
        {
            int nsessi = numt(posdp(i)); /* number of treatments for type i */
            int Mu = nsessi; /* each person initialized with their own cluster */
            /* cluster identifiers */
            slabs = seq_len(Mu);
	    su = as<ucolvec>(slabs);
            su -= 1; su = shuffle(su);
            ss(i,0).set_size(Mu);
            ss(i,0) = su;
            /* cluster locations */
            ustars(i,0) = randn<colvec>(Mu)*sqrt(10);
            /* cluster sizes */
            nums(i,0).set_size(Mu); nums(i,0).ones();
            Ms(i) = Mu;
            /* return objects */
            Sulist(i,0).set_size(nkeep,nsessi);  
            bigSulist(i,0).set_size(nsessi,2);
        } 
    }
    
    if(numcar > 0) /* return CAR parameters for each CAR structure */
    {
        /* set sizes for adjacency matrix (arma) structures */
	ngs = as<icolvec>(numgs);
        if( numcar != nOmega )
                throw std::range_error("CAR entries not equal to number of Omega inputs");
        omegamats.set_size(numcar,1); 
        omegapluses.set_size(numcar,1); 
        for(i = 0; i < numcar; i++)
        {          
                omegamats(i,0) = as<mat>(Omega[i]);
                omegapluses(i,0) = as<rowvec>(omegaplus[i]);
        }
    }

    if(numigp > 0)
    {
        Etalist.set_size(numigp,1);
        etas.set_size(numigp,1);
        tauetas.set_size(numigp); tauetas.ones();
        mmats.set_size(numigp,1);
        for(i = 0; i < numigp; i++)
        {
            /* sampling objects */
            mmats(i,0) = as<mat>(Mmat[i]);
            int ngi = mmats(i,0).n_cols;
            etas(i,0) = randn<colvec>(ngi);
            /* return objects */
            Etalist(i,0).set_size(nkeep,ngi);
        }
    }
    
    // Initialize other (non-session) parameter values
    /* subjects cluster capture variables */
    double taue = 1;
    rowvec taub(nr); taub.ones();
    IntegerVector sv = seq_len(np);
    icolvec s(sv.begin(), np, false);
    s -= 1; s = shuffle(s);
    int M = np; /* each person initialized with their own cluster */
    icolvec num(M); num.ones();
    mat    bstarmat = randn<mat>(M,nr)*sqrt(1/taub(0));
    mat    bmat(np,nr); bmat.zeros();
    colvec zb(nc); zb.zeros();
    double conc = 1; /* DP concentration parameter */
    /* initialize constant, c, used to compose ytilde */
    colvec cmm(nc); cmm.zeros();
    for(int p = 0; p < (nty-1); p++) /* leave out the last treatment - will add back in function 'cu' */
    {
        cmm += wmats(p,0)*us(p,0);
    }
    /* remaining parameters */
    colvec beta = randn<colvec>(nf)*sqrt(0.1);
    colvec resid(nc); resid.zeros(); /* regression residual */
    double alpha = rnorm( 1, 0, sqrt(10) )[0];
    mat pbmat(nr,nr); pbmat.eye();
    /* fit assessment - related measures */
    double deviance = 0;  ucolvec ordscore; ordscore.zeros();
    mat phat(np,np); phat.zeros(); /* empirical co-clustering probability matrix */
    rowvec devmarg(nc); colvec devres(4); devres.zeros();
    rowvec logcpo(nc); logcpo.zeros(); double lpml;
    
    /* capture samples to return */
    int oo = 0, kk;

    // set hyperparameter values
    double a6, b6;
    a6 = ac; b6 = bc; /* conc */
    
    // Conduct posterior samples and store results
    int k,l;
    for(k = 0; k < niter; k++)
    {
	//if( (k % 1000) == 0 ) Rcout << "Interation: " << k << endl;
        clusterbstep(xmat, zmat, wmats, us, zsplit, ytilsplit, pbmat, y, 
            beta, alpha, taue, taub, persons, bstarmat, s, num, M, conc,
            np, nr);
        concstep(conc, M, np, a6, b6);
        bstarstep(zsplit, ytilsplit, taue, pbmat, s, num, bstarmat,
                zb, bmat, np, nr);
        betalsmultstep(xmat, wmats, us, y, beta, alpha, taue, zb, nf);
        umultsteps(us, ustars, nums, ss, concs, Ms, etas, tauus, tauetas, omegamats, 
             omegapluses, wmats, mmats, ngs, typet, xmat, y, beta, zb, alpha, taue, cmm,
             ustrength);
        alphalsmultstep(xmat, wmats, us, beta, zb, y, resid, alpha, taue, nc);
        taumultstep(bstarmat, resid, taub, taue, M, nc);
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
                
                for(int i = 0; i < nty; i++)
                {
                    U(i,0).row(oo) = trans( us(i,0) );
                }
                
                if(numdp > 0)
                {
                    numMu.row(oo) = Ms;
                    Concu.row(oo) = concs;
                    for(int i = 0; i < numdp; i++)
                    {
                        Sulist(i,0).row(oo) = trans( ss(i,0) );
                    }
                } 
                
                if(numigp > 0)
                {
                    for(int i = 0; i < numigp; i++)
                    {
                        Etalist(i,0).row(oo) = trans( etas(i,0) );
                    }
                }
                Tauu.row(oo) = tauus;
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
    
    // compute lease squares cluster for clustering of sessions
    if(numdp > 0)
    {
        ucolvec ordscoreu;
        for(int i = 0; i < numdp; i++)
        {
            lsqclusteru(Sulist(i,0), ordscoreu, bigSulist(i,0));
        }
    }
    
    // compute least squares cluster for clustering of subjects
    lsqcluster(S, Num, ordscore, phat, bigS);
    optPartition(0,0) = phat;
    optPartition(1,0) = conv_to<mat>::from(ordscore);
    optPartition(2,0) = conv_to<mat>::from(S);
    // DIC
    dic3comp(Deviance, Devmarg, devres); /* devres = c(dic,dbar,dhat,pd) */
    cpo(Devmarg, logcpo, lpml);

    // Return results
    return Rcpp::List::create(Rcpp::Named("Deviance") = Deviance,
                                  //Rcpp::Named("Devmarg") = Devmarg,
                                  Rcpp::Named("devres") = devres,
                                  //Rcpp::Named("logcpo") = logcpo,
                                  Rcpp::Named("lpml") = lpml,
				  Rcpp::Named("Beta") = Beta,
				  Rcpp::Named("Alpha") = Alpha,
                                  Rcpp::Named("B") = B,
                                  //Rcpp::Named("C") = Conc,
                                  Rcpp::Named("U") = U,
                                  Rcpp::Named("Tauu") = Tauu,
                                  Rcpp::Named("Su") = Sulist,
                                  Rcpp::Named("Mu") = numMu,
                                  Rcpp::Named("bigSu") = bigSulist,
                                  Rcpp::Named("Cu") = Concu,
                                  Rcpp::Named("Taue") = Taue,
                                  Rcpp::Named("Taub") = Taub,
                                  Rcpp::Named("Residuals") = Resid,
                                  Rcpp::Named("M") = numM,
                                  //Rcpp::Named("S") = S,
                                  Rcpp::Named("Num") = Num,
                                  Rcpp::Named("optPartition") = optPartition,
                                  Rcpp::Named("bigSmin") = bigS
				  );
END_RCPP
} /* end MCMC function returning SEXP */

SEXP umultsteps(field<colvec>& us, field<colvec>& ustars, field<icolvec>& nums,
           field<ucolvec>& ss, rowvec& concs, irowvec& Ms, field<colvec>& etas, 
           rowvec& tauus, rowvec& tauetas, const field<mat>& omegamats, 
           const field<rowvec>& omegapluses, const field<mat>& wmats, 
           const field<mat>& mmats, const icolvec& ngs, const icolvec& typet, 
           const mat& xmat, const colvec& y,  const colvec& beta, const colvec& zb, 
           double alpha, double taue, colvec& cmm, double ustrength)
    {
        BEGIN_RCPP
        int k, sel, countcar = 0, countind = 0, countigrp = 0, countdp = 0;
        int ns, ng;
        double a, b;
        mat qmat;
        int nty = typet.n_elem;
        for( k = 0; k < nty; k++ )
        {
            if(typet(k) == 1) /* mmcar */
            {
                sel = countcar;
                ucarstep(xmat, beta, zb, y, wmats, us, omegamats(sel,0), 
                        omegapluses(sel,0), alpha, taue, tauus(k), k, cmm);
                qmat = -omegamats(sel,0); qmat.diag() = omegapluses(sel,0);
		if(ngs.n_elem > 1)
		{
                	ng = ngs(sel); 
		}else{
			ng = 0;
		}
		ns = wmats(k,0).n_cols;
                a = 0.5*(double(ns) - double(ng)) + ustrength;
                b = 0.5*( as_scalar(trans(us(k,0))*qmat*us(k,0)) ) + ustrength;
                tauus(k) = rgamma(1, a, (1/b))[0];
                countcar += 1;
            }else{
                if(typet(k) == 2) /* ind */
                {
                    sel = countind;
                    uistep(xmat, beta, zb, y, wmats, us, alpha, taue, tauus(k), k, cmm);
                    ns = wmats(k,0).n_cols;
                    a = 0.5*double(ns) + ustrength;
                    b = 0.5* dot( us(k,0) , us(k,0) ) + ustrength;
                    tauus(k) = rgamma(1, a, (1/b))[0];
                    countind += 1;
                }else{
                    if(typet(k) == 3) /* igp */
                    {
                        sel = countigrp;
                        ns = mmats(sel,0).n_rows;
                        ng = mmats(sel,0).n_cols;
                        uigrpstep(xmat, mmats(sel,0), beta, zb, y, etas(sel,0), 
                                wmats, us, alpha, taue, tauus(k), k, cmm);
                        etagrpstep(mmats(sel,0), etas(sel,0), us(k,0),
                                tauus(k,0), tauetas(sel));
                        /* tauu */
                        a = 0.5*double(ns)+ ustrength;
                        colvec uresid = us(k,0) - mmats(sel,0)*etas(sel,0);
                        b = 0.5* dot( uresid , uresid ) + ustrength;
                        tauus(k) = rgamma(1, a, (1/b))[0];
                        /* taueta */
                        a = 0.5*double(ng) + 1;
                        b = 0.5* dot( etas(sel,0) , etas(sel,0) ) + 1;
                        tauetas(sel) = rgamma(1, a, (1/b))[0];
                        countigrp += 1;
                    }else{ /* typet(k) is dp */
                        sel = countdp;
                        int nc = xmat.n_rows;
                        ns = wmats(k,0).n_cols;
                        colvec ytilde(nc); ytilde.zeros();
                        clusterustep(y, zb, xmat,  beta, alpha, taue, tauus(k),
                                wmats, us, ustars(sel,0), ss(sel,0), nums(sel,0), 
                                Ms(sel), concs(sel), ytilde, k, cmm);
                        double a6 = 1, b6 = 1;
                        concstep(concs(sel), Ms(sel), ns, a6, b6);
                        ustarstep(ytilde, wmats(k,0), ustars(sel,0), us, 
                                ss(sel,0), taue, tauus(k), k);
                        a = 0.5*double(Ms(sel)) + 1;
                        b = 0.5* dot( ustars(sel,0) , ustars(sel,0) ) + 1;
                        tauus(k) = rgamma(1, a, (1/b))[0];
                        countdp += 1;
                    }
                }
            }
        }
        
        END_RCPP
    }

    SEXP clusterustep(const colvec& y, const colvec& zb, const mat& xmat,
            const colvec& beta, double alpha, double taue, double tauu,
            const field<mat>& wmats, const field<colvec>& us, colvec& ustar,
            ucolvec& s, icolvec& num, int& M, double& conc, colvec& ytilde,
	    int treat, colvec& cmm)
    {
        BEGIN_RCPP
        // compute number of cases for each person and store results in
        int ns = wmats(treat,0).n_cols;
        double logq0, q0, sweights;
        double ej, hj, phij, ustarj;
        int j, l, m;
        colvec cj, wj, ytildej, cstarj, usj;
        mat wlessj;
        
        /* response less nuisance parameters */
        colvec c; cu(zb, xmat, beta, alpha, wmats, us, cmm, c, treat);     
        /* return data object for location sampling */
	ytilde = y - c;  /* return for posterior sampling of locations, ustar, in a joint update */
            
        for(j = 0; j < ns; j++)
        {
            /* remove all but focus session parameter from y */
            wlessj = wmats(treat,0); usj = us(treat,0);
            wlessj.shed_col(j);
            usj.shed_row(j);
            cj = c + wlessj*usj;
            ytildej = y - cj;

            // sample cluster assignment indicators, s(1),...,s(np)
            if(num(s(j)) == 1) /* remove singleton cluster */
            {
                ustar.shed_row(s(j));
                num.shed_row(s(j));
                M -= 1; /* decrement cluster count */

                //decrement cluster tracking values by 1 for tossed cluster
                uvec ord = sort_index(s); /* of size np, same as s */
                int k; int thresh = ns;
                for (k = 0; k < ns; k++)
                {
                    if( s(ord(k)) > s(j) )
                    {
                        thresh = k;
                        break;
                    }
                } /* end loop k to find thresh*/
                // decrement s > s(j)
                for (m = thresh; m < ns; m++)
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
            NumericVector zro(1); zro(0) = 0;
            logq0 = loglike(ytildej,taue);
            logq0 += dnorm(zro, 0.0, sqrt(1/tauu), true)[0]; /* true = log */
            wj = wmats(treat,0).col(j);
            ej = taue*dot(wj,ytildej);
            phij = taue*dot(wj,wj) + tauu; 
            hj = ej/phij;
            logq0 -= dnorm(zro, hj, sqrt(1/phij), true)[0]; /* using symmetry to switch hj and zro */
            q0 = exp(logq0);

            // construct posterior sampling weights to sample s(j)
            colvec weights(M+1); weights.zeros();
            for (l = 0; l < M; l++) /* cycle through all clusters for s(j) */
            {
                s(j) = l; /* will compute likelihoods for every cluster */
                ustarj = ustar.row(s(j));
                cstarj = cj + wj*ustarj;
                ytildej = y - cstarj;

                weights(l) = exp(loglike(ytildej,taue));
                weights(l) *= num(s(j))/(double(ns) - 1 + conc);
            }
            weights(M) = conc/(double(ns)-1+conc)*q0;

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
            if(s(j) == (unsigned)M)
            {
                ustarj = rnorm(1,hj,sqrt(1/phij))[0];
                ustar.insert_rows(M,1);
                num.insert_rows(M,1);
                ustar.row(M) = ustarj;
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

    SEXP ustarstep(const colvec& ytilde, const mat& wtreat,
	    colvec& ustar, field<colvec>& us, const ucolvec& s, 
	    double taue, double tauu, double treat)
    {
        BEGIN_RCPP
        // sample ustar, jointly
        int M = ustar.n_elem;
 	int ns = wtreat.n_cols;
	int nc = ytilde.n_elem;

        mat C(ns,M); C.zeros();
	for(int i = 0; i < ns; i++)
        {
		C(i,s(i)) = 1;
	}
        mat pmat = tauu*eye<mat>(M,M);
	colvec cons(nc); cons.zeros();
        rmvnchol(wtreat*C, pmat, ytilde, cons, ustar, M, taue);
        us(treat,0) = ustar.elem(s);

        END_RCPP
    } /* end function bstep for sampling bmat and zb */



    SEXP lsqclusteru(const umat& S, ucolvec& ordscore, umat& bigS)
    {
        BEGIN_RCPP
        int n = S.n_rows;  int ns = S.n_cols;
        mat delta; delta.eye(ns,ns);
        mat phat(ns,ns); phat.zeros();
        mat diffsq(ns,ns); diffsq.zeros();
        urowvec s(ns); colvec score(n);
	ucolvec ord(ns); ord.zeros();
	unsigned long minscoreloc;
        int i, j, k;
        // create phat
        for(i = 0; i < n; i++)
        {
	    delta.eye();
            s = S.row(i);
            for(j = 0; j < ns; j++)
            {
                for(k = j + 1; k < ns; k++)
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
            for(j = 0; j < ns; j++)
            {
                for(k = j + 1; k < ns; k++)
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
        urowvec soptim = S.row(minscoreloc);
        bigS.set_size(ns,2);
        for(int i = 0; i < ns; i++)
        {
            bigS(i,0) = i + 1; /* adding 1 shifts the subject index by 1 for R */
            bigS(i,1) = soptim(i);
                    
        }

        END_RCPP
    } /* end function lsqclusteru */

    SEXP ucarstep(const mat& xmat, const colvec& beta, const colvec& zb,
            const colvec& y, const field<mat>& wmats, field<colvec>& us,
	    const mat& omega, const rowvec& omegaplus,
            double alpha, double taue, double tauu, int treat, colvec& cmm)
    {
        BEGIN_RCPP
        // set up structures that will set column s of wcase = 0 and
        // entry s for row s of omega = 0 for sampling u[s]
	int ns = omega.n_rows;
        int nc = xmat.n_rows;
        colvec ws(nc), ytilde(nc);
        mat wmat = wmats(treat,0);
	colvec u = us(treat,0);
        double es, hs, phis;        

	// build portion of offset constant, c, not dependent on u = us(treat,0)
	colvec cwithalls; cu(zb, xmat, beta, alpha, wmats, us, cmm, cwithalls, treat);
	// S x 1 objects
	cwithalls 		+= wmat*u; /* cu adds in previously sampled MM term and removes currently sampled.  Add it back in */
	colvec omega_uall	= omega*u;
        // loop over s to sample u[s]|u[-s],..,y
        for(int s = 0; s < ns; s++)
        {
            // Set entry s for data to 0
            ws = wmat.col(s); /* w.s = W.mat[,s] */
            // remove entry W.zeros[,s]*u[s] from cwithalls
	    // removing S x 1 objects based on u(s)
	    cwithalls 	-= ws*u(s);
	    omega_uall	-= omega.col(s)*u(s);
            // sample u[s]
            // construct posterior mean, hs, and precision, phis
            ytilde 	= y - cwithalls;
            es 		= taue*dot(ytilde,ws) + tauu*omega_uall(s);
            phis 	= taue*dot(ws,ws) + tauu*omegaplus(s);
            hs 		= es*(1/phis);
            u(s) 	= rnorm( 1, hs, sqrt(1/phis) )[0];
	    // puts back session s contribution from newly sampled value for u(s)
            cwithalls 	+= ws*u(s);
	    omega_uall	+= omega.col(s)*u(s);
        } /* end loop over s for sampling u */
        // post-process u to subtract out group means since CAR prior
        // identified only up to difference in means
        u -= mean(u);
        us(treat,0) = u;

        END_RCPP
    } /* end function ustep to sample u[1,...,s] */
    
    SEXP uistep(const mat& xmat, const colvec& beta, const colvec& zb,
            const colvec& y, const field<mat>& wmats, field<colvec>& us,  
            double alpha, double taue, double tauu, int treat, colvec& cmm)
    {
        BEGIN_RCPP
        // Sample posterior of u (ns x 1) from posterior Gaussian
        int ns = wmats(treat,0).n_cols;
        colvec u(ns);
        // compute prior precision matrix
        mat pu(ns,ns); pu.eye(); pu *= tauu;
        colvec c; cu(zb, xmat, beta, alpha, wmats, us, cmm, c, treat); /* nc x 1 */
        rmvnchol(wmats(treat,0), pu, y, c, u, ns, taue);
        us(treat,0) = u;

        END_RCPP
    } /* end function to sample nf x 1 fixed effects, beta */
    
    SEXP uigrpstep(const mat& xmat, const mat& mmat, const colvec& beta, 
            const colvec& zb, const colvec& y, const colvec& eta, 
            const field<mat>& wmats, field<colvec>& us, double alpha, double taue,
            double tauu, int treat, colvec& cmm)
    {
        BEGIN_RCPP
        /* memo: must use eta, rather than etas since treat references to nty, not numigrp */
        // Sample posterior of u (ns x 1) from posterior Gaussian
        int ns = wmats(treat,0).n_cols;
        colvec u(ns);
        // compute prior precision matrix
        mat pu(ns,ns); pu.eye(); pu *= tauu;
        colvec c; cu(zb, xmat, beta, alpha, wmats, us, cmm, c, treat); /* nc x 1 */
        colvec m = mmat*eta; /* ns x 1 */
        rmvnmean(wmats(treat,0), pu, y, c, m, u, ns, taue);
        us(treat,0) = u;

        END_RCPP
    } /* end function to sample nf x 1 fixed effects, beta */
    
    SEXP etagrpstep(const mat& mmat, colvec& eta, const colvec& u,
            double tauu, double taueta)
    {
        BEGIN_RCPP
        // Sample posterior of u (ns x 1) from posterior Gaussian
        // compute prior precision matrix
        int ns = u.n_elem, ng = eta.n_elem;
        mat pu(ng,ng); pu.eye(); pu *= taueta;
        colvec cons(ns); cons.zeros();
        rmvnchol(mmat, pu, u, cons, eta, ng, tauu);
        // enforce constraint for mean(eta) = 0
        eta -= mean(eta);

        END_RCPP
    } /* end function to sample nf x 1 fixed effects, beta */

    SEXP cu(const colvec& zb, const mat& xmat,
            const colvec& beta, double alpha, 
            const field<mat>& wmats, const field<colvec>& us,
	    colvec& cmm, colvec& c, int treat)
    {
    	BEGIN_RCPP
	int nty = wmats.n_rows;
	c = alpha + xmat*beta + zb;
        /* method limits number of multiplications per cmm iteration to 2 */
        if(treat == 0) // first treatment, k = 0
	{
		cmm += wmats((nty-1),0)*us((nty-1),0);
	}else{
		cmm += wmats((treat-1),0)*us((treat-1),0);
	}
	cmm -= wmats(treat,0)*us(treat,0);
        c += cmm;
    	END_RCPP
    }

    SEXP clusterbstep(const mat& xmat, const mat& zmat, const field<mat>& wmats,
            const field<colvec>& us, field<mat>& zsplit,
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
        int nty = wmats.n_rows;
        field<mat> wjs(nty,1);
        mat zj, xj, yj, phib(nr,nr), phibinv(nr,nr);
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
            for(l = 0; l < nty; l++)
            {
                wjs(l,0) = wmats(l,0).rows(startrow,endrow);
                cj += wjs(l,0)*us(l,0);
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
                logq0 += dnorm(zro, 0.0, sqrt(1/taub(i)), true)[0]; /* true = log */
            }
            eb = taue*trans(zj)*ytildej;
            phib = taue*trans(zj)*zj + pbmat; phibinv = inv(phib);
            hb = phibinv*eb;
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

    SEXP betalsmultstep(const mat& xmat, const field<mat>& wmats, 
                const field<colvec>& us, const colvec& y,
                colvec& beta, double alpha, double taue,
                const colvec& zb, int nf)
    {
        BEGIN_RCPP
        // Sample posterior of beta (nf x 1) from posterior Gaussian
        int nty = wmats.n_rows;
        colvec c = alpha + zb; /* nc x 1 */
        for(int i = 0; i < nty; i++)
        {
            c += wmats(i,0)*us(i,0);
        }
        rmvnlschol(xmat, y, c, beta, nf, taue);
        END_RCPP
    } /* end function to sample nf x 1 fixed effects, beta */
    
    SEXP alphalsmultstep(const mat& xmat, const field<mat>& wmats, 
                const field<colvec>& us, const colvec& beta,
                const colvec& zb, const colvec& y, 
                colvec& resid, double& alpha, double taue,
                long nc)
    {
        BEGIN_RCPP
        int nty = wmats.n_rows;
        colvec c = xmat*beta + zb; /* nc x 1 */
        for(int i = 0; i < nty; i++)
        {
            c += wmats(i,0)*us(i,0);
        }
        resid = y - c;
        double ea = taue*sum(resid);
        double phia = taue*double(nc);
        double ha = ea*(1/phia);
        alpha = rnorm( 1, ha, sqrt(1/phia) )[0];
        resid -= alpha;
        END_RCPP
    } /* end function to sample intercept, alpha */
    
    SEXP taumultstep(const mat& bstarmat, const colvec& resid, 
            rowvec& taub, double& taue, int M, int nc)
    {
        BEGIN_RCPP
        double a, b;
        /* taub */
        int nr = bstarmat.n_cols;
        a = 0.5*double(M) + 1;
        int i;
        for(i = 0; i < nr; i++)
        {
            b =  0.5*dot( bstarmat.col(i), bstarmat.col(i) ) + 1;
            taub(i) = rgamma(1, a, (1/b))[0];
        }
        /* taue */
        a = 0.5*double(nc) + 1;
        b = 0.5*dot( resid, resid ) + 1;
        taue = rgamma(1, a, (1/b))[0];
        END_RCPP
    } /* end function taustep to sample precision parameters */
