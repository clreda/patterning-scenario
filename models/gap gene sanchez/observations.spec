// Trajectory

#Trajectory[0].GRNs{1} |= $InitialGRNs1;
#Trajectory[0].Phenotypes{1} |= $InitialPhenotypes1;
#Trajectory[0].GRNs{2} |= $InitialGRNs1;
#Trajectory[0].Phenotypes{2} |= $InitialPhenotypes1;
#Trajectory[0].GRNs{3} |= $InitialGRNs1;
#Trajectory[0].Phenotypes{3} |= $InitialPhenotypes2;
#Trajectory[0].GRNs{4} |= $InitialGRNs1;
#Trajectory[0].Phenotypes{4} |= $InitialPhenotypes2;
#Trajectory[5].Phenotypes{1} |= $FinalPhenotypes1;
#Trajectory[5].Phenotypes{2} |= $FinalPhenotypes2;
#Trajectory[5].Phenotypes{3} |= $FinalPhenotypes3;
#Trajectory[5].Phenotypes{4} |= $FinalPhenotypes4;
fixpoint(#Trajectory[5].Phenotypes{1});
fixpoint(#Trajectory[5].Phenotypes{2});
fixpoint(#Trajectory[5].Phenotypes{3});
fixpoint(#Trajectory[5].Phenotypes{4});

// Gene Regulatory Networks

$InitialGRNs1 := {
	$Gt := {
		$formula := {
			( ( ( Hb-zyg=0 ) and ( Kr=0 ) and Bcd=1 ) or ( ( Hb-zyg=0 ) and ( Kr=0 ) and Cad=2 ) or ( ( Hb-zyg=0 ) and Bcd=1 and Cad=2 ) or ( ( Kr=0 ) and Bcd=1 and Cad=2 ) );
			$output := 1;
		};
	}
	$Hb-zyg := {
		$formula := {
			( ( Hb-zyg=2 and ( Kr=0 or Kr=1 ) ) or ( Hb-zyg=3 and ( Kr=0 or Kr=1 ) ) or ( Hb-mat=1 and ( Kr=0 or Kr=1 ) ) or ( Bcd=1 and ( Kr=0 or Kr=1 ) ) or ( Bcd=2 and ( Kr=0 or Kr=1 ) ) or ( Bcd=3 and ( Kr=0 or Kr=1 ) ) );
			$output := 1;
		};
		$formula := {
			( ( Hb-zyg=1 and ( Kr=0 or Kr=1 ) ) or ( Hb-zyg=3 and ( Kr=0 or Kr=1 ) ) or ( Hb-mat=1 and ( Kr=0 or Kr=1 ) ) or ( Bcd=1 and ( Kr=0 or Kr=1 ) ) or ( Bcd=2 and ( Kr=0 or Kr=1 ) ) or ( Bcd=3 and ( Kr=0 or Kr=1 ) ) );
			$output := 2;
		};
		$formula := {
			( ( Hb-zyg=1 and Hb-zyg=2 and Bcd=1 and Bcd=3 ) );
			$output := 3;
		};
	}
	$Hb-mat := {
		$formula := {
			( ( Hb-mat=2 ) or ( Hb-mat=3 ) );
			$output := 1;
		};
		$formula := {
			( ( Hb-mat=1 ) or ( Hb-mat=3 ) );
			$output := 2;
		};
		$formula := {
			( ( Hb-mat=1 ) or ( Hb-mat=2 ) );
			$output := 3;
		};
	}
	$Kr := {
		$formula := {
			( ( Hb-zyg=1 and Kr=2 and Bcd=1 ) );
			$output := 1;
		};
		$formula := {
			( ( Hb-zyg=1 and ( Gt=0 ) ) or ( Kr=1 and ( Gt=0 ) ) or ( Bcd=1 and ( Gt=0 ) ) or ( Hb-zyg=1 and Kr=1 and Bcd=1 ) );
			$output := 2;
		};
	}
	$Kni := {
		$formula := {
			( ( Bcd=1 and ( Gt=0 ) and ( Hb-zyg=0 or Hb-zyg=1 ) ) or ( Cad=1 and ( Gt=0 ) and ( Hb-zyg=0 or Hb-zyg=1 ) ) );
			$output := 1;
		};
	}
	$Bcd := {
		$formula := {
			( ( Bcd=2 and Bcd=3 ) );
			$output := 1;
		};
		$formula := {
			( ( Bcd=1 ) or ( Bcd=3 ) );
			$output := 2;
		};
		$formula := {
			( ( Bcd=1 ) or ( Bcd=2 ) );
			$output := 3;
		};
	}
	$Cad := {
		$formula := {
			( ( Cad=2 ) );
			$output := 1;
		};
		$formula := {
			( ( Cad=1 ) );
			$output := 2;
		};
	}
};

// Initial phenotypes

$InitialPhenotypes1 := {
	Gt = 0 and
	Hb-zyg = 2 and
	Kr = 0 and
	Kni = 0
};

$InitialPhenotypes2 := {
	Gt = 0 and
	Hb-zyg = 0 and
	Kr = 0 and
	Kni = 0
};

// Intermediate phenotypes

$IntermediatePhenotypes11 := {
	Gt = 1 and
	Hb-zyg = 2 and
	Kr = 0 and
	Kni = 0
};

$IntermediatePhenotypes12 := {
	Gt = 0 and
	Hb-zyg = 3 and
	Kr = 0 and
	Kni = 0
};

$IntermediatePhenotypes21 := {
	Gt = 0 and
	Hb-zyg = 2 and
	Kr = 1 and
	Kni = 0
};

$IntermediatePhenotypes31 := {
	Gt = 0 and
	Hb-zyg = 1 and
	Kr = 0 and
	Kni = 0
};

$IntermediatePhenotypes32 := {
	Gt = 0 and
	Hb-zyg = 0 and
	Kr = 0 and
	Kni = 1
};

$IntermediatePhenotypes33 := {
	Gt = 0 and
	Hb-zyg = 1 and
	Kr = 0 and
	Kni = 1
};

$IntermediatePhenotypes34 := {
	Gt = 0 and
	Hb-zyg = 1 and
	Kr = 1 and
	Kni = 0
};

$IntermediatePhenotypes41 := {
	Gt = 1 and
	Hb-zyg = 0 and
	Kr = 0 and
	Kni = 1
};

$IntermediatePhenotypes42 := {
	Gt = 0 and
	Hb-zyg = 0 and
	Kr = 0 and
	Kni = 1
};

// Final phenotypes

$FinalPhenotypes1 := {
	Gt = 1 and
	Hb-zyg = 3 and
	Kr = 0 and
	Kni = 0
};

$FinalPhenotypes2 := {
	Gt = 0 and
	Hb-zyg = 2 and
	Kr = 2 and
	Kni = 0
};

$FinalPhenotypes3 := {
	Gt = 0 and
	Hb-zyg = 1 and
	Kr = 1 and
	Kni = 1
};

$FinalPhenotypes4 := {
	Gt = 1 and
	Hb-zyg = 0 and
	Kr = 0 and
	Kni = 0
};
