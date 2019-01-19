
carstm_hyperparameters = function( reference_sd, alpha=0.5, reference_mean=0 ) {
  # some generic PC priors, scaled by sd of data
  # pc.prior to median .. minimally info. scale

  hyper = list(

    iid = list(
      prec = list(
        prior = "pc.prec",  # exponential decay
        param = c(reference_sd, alpha)
      )
    ),

    # means informative, sd marginally diffuse
    # see: inla.set.control.fixed.default() for defaults
    fixed = list(
        mean.intercept = reference_mean,
        prec.intercept = 1e-3,
        mean=0,
        prec=1e-2
    ),


    # param=c(u, alpha); u=sigma; alpha=prob;
    # see inla.doc("pc.rw2") inla.doc("pc.prec")  ..prior sd attributable to rw2
    rw2 = list(
      prec = list(
        prior = "pc.prec",  # exponential decay
        param = c(reference_sd, alpha)
      )
    ),

    # see inla.doc("ar1") ; theta0, theta1 are expected
    # param=c(u, alpha); u=sigma; alpha=prob;
    # see inla.doc("pc.prec")  ..prior sd attributable to autocor rho
    # param=c(u, alpha); rho = 0.5; u=sqrt(1-rho); alpha=prob; see inla.doc("pc.cor1")
    ar1 = list(
      prec = list(
        prior = "pc.prec",  # exponential decay
        param = c(reference_sd, alpha)
      ),
      rho = list(
        prior = "pc.cor0", # inla.doc("pc.cor0") ..base model: rho = 0  --- expoential; will tend to 0 unless there is info
        param = c(sqrt(1-0.5), 0.1)  # rho=0.5; u=sqrt(1-rho)  ... 100-10% of probablity weight to rho 0.5 or less .. forces smooth and only goes high if really high
      )
    ),

    # param=c(u, alpha); u=phi (proportion spatial); alpha=prob
    bym2 = list(
      prec = list(
        prior = "pc.prec",
        param = c(reference_sd, alpha)
      ),
      phi = list(
        prior="pc",  # see bottom of inla.doc("bym2")
        param=c(0.5, 0.5) # c(phi=0.5, alpha=0.5)
      )
    )
  )

  return(hyper)
}
