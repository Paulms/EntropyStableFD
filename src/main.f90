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
  USE test2
  USE test3
  USE test4
  IMPLICIT NONE

  ! Run Test(run_reference, run_error)
  ! CALL burger_runexample(1, .true., .true.)
  ! CALL test2_run(.false., .true.)
  ! CALL test3_run(.false., .false.)
  CALL test4_run()

END PROGRAM EntropyStableFD