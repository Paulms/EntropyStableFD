PROGRAM EntropyStableFD
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Entropy Stable Schemes for Degenerate Equations      !!
  !!                                                      !!
  !!                                                      !!
  !! Author: Paul Mendez Silva                            !!
  !! e-mail: paul.mendez@udec.cl                          !!
  !! Date: 10/Marzo/2016                                  !!
  !!                                                      !!
  !! Version: 0.1                                         !!  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE example_burger
  IMPLICIT NONE
  CALL burger_runexample(1)

END PROGRAM EntropyStableFD