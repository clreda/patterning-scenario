// Trajectory

#Trajectory[0].GRNs{1} |= $InitialGRNs1;
#Trajectory[0].GRNs{2} |= $InitialGRNs2;
#Trajectory[0].GRNs{3} |= $InitialGRNs1;
#Trajectory[0].GRNs{4} |= $InitialGRNs1;
#Trajectory[0].Phenotypes{1} |= $InitialPhenotypes;
#Trajectory[0].Phenotypes{2} |= $InitialPhenotypes;
#Trajectory[0].Phenotypes{3} |= $InitialPhenotypes;
#Trajectory[0].Phenotypes{4} |= $InitialPhenotypes;
#Trajectory[10].Phenotypes{1} |= $FinalPhenotypes1;
#Trajectory[10].Phenotypes{2} |= $FinalPhenotypes2;
#Trajectory[10].Phenotypes{3} |= $FinalPhenotypes2;
fixpoint(#Trajectory[10].Phenotypes{1});
fixpoint(#Trajectory[10].Phenotypes{2});
fixpoint(#Trajectory[10].Phenotypes{3});

// Gene Regulatory Network

$InitialGRNs1 := {
	$M1 := {
		$formula := {
			( M1=1 and M2=0 );
			$output := 1;
		};
		$formula := {
			( ( M1=1 or M1=2 ) and M2=0 );
			$output := 2;
		};
	}
	$M2 := {
		$formula := {
			( M2=1 and M1=0 );
			$output := 1;
		};
	}
	$A := {
		$formula := {
			( ( M1=1 or M2=0 ) and C=0 );
			$output := 1;
		};
	}
	$B := {
		$formula := {
			( A=0 or C=1 );
			$output := 1;
		};
	}
	$C := {
		$formula := {
			( ( M1=1 or M2=0 ) and B=1 );
			$output := 1;
		};
	}
};

$InitialGRNs2 := {
	$M1 := {
		$formula := {
			( M1=1 and M2=0 );
			$output := 1;
		};
		$formula := {
			( ( M1=1 or M1=2 ) and M2=0 );
			$output := 2;
		};
	}
	$M2 := {
		$formula := {
			( M2=1 and M1=0 );
			$output := 1;
		};
	}
	$A := {
		$formula := {
			( ( M1=1 ) and C=0 );
			$output := 1;
		};
	}
	$B := {
		$formula := {
			( A=0 or C=1 );
			$output := 1;
		};
	}
	$C := {
		$formula := {
			( ( M2=0 ) and B=1 );
			$output := 1;
		};
	}
};

$InitialGRNs3 := {
	$M1 := {
		$formula := {
			( M1=1 );
			$output := 1;
		};
		$formula := {
			( ( M1=1 or M1=2 ) and M2=0 );
			$output := 2;
		};
	}
	$M2 := {
		$formula := {
			( M2=1 and M1=0 );
			$output := 1;
		};
	}
	$A := {
		$formula := {
			( ( M1=1 or M2=0 ) and C=0 );
			$output := 1;
		};
	}
	$B := {D
		$formula := {
			( A=0 or C=1 );
			$output := 1;
		};
	}
	$C := {
		$formula := {
			( ( M1=1 or M2=0 ) and B=1 );
			$output := 1;
		};
	}
};
// Initial phenotypes

$InitialPhenotypes := {
	A = 0 and
	B = 0 and
	C = 0
};

// Final phenotypes

$FinalPhenotypes1 := {
	A = 0 and
	B = 1 and
	C = 0
};

$FinalPhenotypes2 := {
	A = 0 and
	B = 1 and
	C = 1
};
