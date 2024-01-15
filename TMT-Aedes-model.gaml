/**
* Name: TMT-Aedes-model
* Author: Samuel Beach 
*/


model TMT

global {
	
//	Experimental variables
	float DDM_adj <- 0.01;
	float tx_competitiveness <- 0.031;
	float release_ratio <- 25.0;
	int release_period <- 3;
	int release_delay <- 365;
	float TMT_lethality <- 1.00;
	float TMT_reduced_lethality <- 0.55;
	
//	Global time variables
	int day <- 24;
	int week <- day * 7;
	int year <- day * 365;
	int max_days <- year * 2;
	
//	Initialisation parameters
	int nb_init <- 2500;
	bool is_batch <- false;
	
//	Immature stage-specific parameters
	float trans_surv <- 0.7;
	float imm_surv <- 0.9995;
	int nb_larvae update: length(larvae) + length(larvae_L3);
	float L1L2_surv function: nb_larvae > 0? imm_surv / (nb_larvae^DDM_adj): 0.9995;
	
//	Mosquito-specific parameters
	float biting_success <- 0.8;
	float adult_surv <- 0.98;
	float mating_prob <- 0.95;
	float fem_delta <- 0.005;
	float male_delta <- 0.012;
	
//	Tracking female parameters
	int mate1 <- 0;
	int mate2 <- 0;
	int mate3 <- 0;
	int mateplus <- 0;
	int bite_count;
	int bite_plus_count;
	int p50 <- 0;
	int p95 <- 0;
	int GC_plus_len <- 0;
	int GC_plus_init;
	
//	transgenic parameters
	bool exp_fsRIDL <- false;
	bool exp_TMT_reduced <- false;
	bool exp_TMT <- false;
	bool exp_SIT <- false;
	
//	World time, and daytime progression
	float step <- 1 #hours;
	int cur_hour update: cycle mod day;
	int nb_days update: int(cycle / day);
	bool new_day function: cur_hour = 0 and cycle != 0? true: false;
	bool night function: cur_hour <= 5 or cur_hour >= 20? true: false;
	bool morning function: cur_hour >= 6 and cur_hour <= 9? true: false;
	bool midday function: cur_hour >= 10 and cur_hour <= 15? true: false;
	bool afternoon function: cur_hour >= 16 and cur_hour <= 19? true: false;
	
//	Release calculations
	int nb_releases <- 0; 
	int sum_Tx <- 0;
	int nb_females_init;
	int nb_moz_release;
	int max_moz_count;
	int daily_moz;
	int total_moz;
	int mean_moz <- 1;
	float surviving_females function: nb_days >= release_delay?
										length(female_mosquito) / nb_females_init:
										0;
	int days_since_release function: nb_days >= release_delay?
										nb_days - release_delay:
										0;
	
	reflex max_count when: (length(male_mosquito) + length(female_mosquito)) > max_moz_count {
		max_moz_count <- (length(male_mosquito) + length(female_mosquito));
	}
	
	reflex mean_count when: new_day and nb_days <= release_delay {
		daily_moz <- length(female_mosquito) + length(male_mosquito);
		total_moz <- total_moz + daily_moz;
		mean_moz <- round(total_moz / nb_days);
	}
	
	reflex release when: new_day and
						 nb_days >= release_delay and 
						 nb_days mod release_period = 0 {
						 	
		nb_moz_release <- round(max_moz_count * release_ratio);
		
		create transgenic_mosquito number: nb_moz_release {
			sum_Tx <- sum_Tx + 1;
		}
		nb_releases <- nb_releases + 1;
	}
	
//	Time to 50% female population reduction
	reflex when: nb_days >= release_delay and p50 = 0 and length(female_mosquito) <= nb_females_init * 0.5 {
		p50 <- cycle - (release_delay*24);
		save [	self.name,
				exp_TMT,
				exp_fsRIDL,
				exp_SIT,
				exp_TMT_reduced,
				p50
		]
		to: "TMT_PR50.csv" format: csv rewrite: false;
	}
	
//	Time to 95% female population reduction
	reflex when: nb_days >= release_delay and p95 = 0 and length(female_mosquito) <= nb_females_init * 0.05 {
		p95 <- cycle - (release_delay*24);
		save [	self.name,
				exp_TMT,
				exp_fsRIDL,
				exp_SIT,
				exp_TMT_reduced,
				p95
		]
		to: "TMT_PR95.csv" format: csv rewrite: false;
	}
	
	reflex step_save when: is_batch and nb_days >= release_delay  {
		save [	self.name,
				self.cycle,
				length(female_mosquito),
				self.nb_females_init,
				GC_plus_len,
				GC_plus_init,
				self.bite_count,
				self.bite_plus_count,
				self.sum_Tx,
				exp_TMT,
				exp_fsRIDL,
				exp_SIT,
				exp_TMT_reduced,
				release_period,
				release_ratio,
				tx_competitiveness
		]
		to: "TMT_step_output.csv" format: csv rewrite: false;
	}
	
	reflex init_female_count when: cycle = (release_delay * 24) {
		nb_females_init <- length(female_mosquito);
		GC_plus_init <- GC_plus_len;
	}

	reflex end when: (nb_days > release_delay) and length(female_mosquito) = 0 and !is_batch {
		do pause;
	}
	
	init {
		create egg number: nb_init {
			age <- round(gauss(4,0.5));
		}
	}
}

species egg  {
	int age;
	int hours_hatch;
	bool fsRIDL_carrier;
	bool male;
	bool female;
	
	init {
		age <- cycle;
		hours_hatch <- round(gauss(4,0.5));
		if (flip(0.5)) {
			female <- true;
		}
		else {
			male <- true;
		}
	}
	
	reflex age {		
		if (!flip(imm_surv)) {
			do die;
		}
		if (cycle >= age + hours_hatch) {
			if (!flip(trans_surv)) {
				do die;
			}
			else {
				create (larvae) {
					self.fsRIDL_carrier <- myself.fsRIDL_carrier;
					self.male <- myself.male;
					self.female <- myself.female;
				}
				do die;
			}
		}
	}
}

species larvae {
	int age;
	int hours_ecdysis;
	bool fsRIDL_carrier;
	bool male;
	bool female;
	float surv;
	
	init {
		age <- cycle;
		hours_ecdysis <- round(gauss(46,6));
		surv <- L1L2_surv;
	}

	reflex age {
		if (L1L2_surv < surv) {
			surv <- L1L2_surv;
		}
		if (!flip(surv)) {
			do die;
		}
		if (cycle >= age + hours_ecdysis) {
			create (larvae_L3) {
				self.fsRIDL_carrier <- myself.fsRIDL_carrier;
				self.male <- myself.male;
				self.female <- myself.female;
			}
			do die;
		}
	}
}

species larvae_L3 {
	int age;
	int hours_pupate;
	bool fsRIDL_carrier;
	bool male;
	bool female;
	
	init {
		age <- cycle;
		hours_pupate <- round(gauss(64,6));
	}

	reflex age {
		if (!flip(imm_surv)) {
			do die;
		}
		if (cycle >= age + hours_pupate) {
			if (!flip(trans_surv)) {
				do die;
			}
			if (female and fsRIDL_carrier) {
				do die;
			}
			else {
				create (pupae) {
					self.fsRIDL_carrier <- myself.fsRIDL_carrier;
					self.male <- myself.male;
					self.female <- myself.female;
				}
				do die;
			}
		}
	}
}

species pupae {
	int age;
	int hours_eclose;
	bool fsRIDL_carrier;
	bool male;
	bool female;
	
	init {
		age <- cycle;
		hours_eclose <- round(gauss(48,6));
	}
	
	reflex age {
		if (!flip(imm_surv)) {
			do die;
		}
		if (cycle >= age + hours_eclose) {
			if (!flip(trans_surv)) {
				do die;
			}
			else if (female) {
				create (female_mosquito) {
					self.fsRIDL_carrier <- myself.fsRIDL_carrier;
				}
			}
			else {
				create (male_mosquito) {
					self.fsRIDL_carrier <- myself.fsRIDL_carrier;
				}
			}
			do die;
		}
	}
}

species mosquito {
	int age <- 0;
	bool SIT <- false;
	bool fsRIDL <- false;
	bool fsRIDL_carrier <- false;
}

species female_mosquito parent: mosquito {
	bool virgin <- true;
	bool gravid <- false;
	bool mated <- false;
	bool fed <- false;
	bool hungry <- false;
	float remating;
	float refeeding <- 0.0006;
	int mated_timer <- 0;
	int num_mates <- 0;
	int num_eggs <- 0;
	int num_bites <- 0;
	int num_GC <- 1;
	bool toxin <- false;
	string death_cause <- "Natural causes";
	int hours_feed;
	int hours_ovi;
	int hours_refeed;
	int hours_since_ovi;
	int hours_since_feed <- 0;
	int egg_timer <- 0;
	int age_hours;
	int ovi_eggs;
	bool has_bitten;
	
	init {
		hours_feed <- round(gauss(36,6));
		hours_ovi <- round(gauss(72,12));
		hours_refeed <- round(gauss(12,4));
		ovi_eggs <- round(gauss(60,10));
		age_hours <- cycle;
	}
	
	reflex polyandry when: mated and (mated_timer > 0 and mated_timer < 25) {
		remating <- (-0.18*ln(mated_timer) + 0.20);
		if (flip(remating)) {
			mated <- false;
		}
		if (mated_timer = 24) {
			remating <- 0.0;
		}
		mated_timer <- mated_timer + 1;
	}
	
	reflex refeeding when: fed {
		if (flip(hours_since_feed * refeeding)) {
			hungry <- true;
		}
	}
	
	reflex age when: new_day{
		age <- age + 1;
		
		float surv <- (-fem_delta * (age)) + (1 + fem_delta);
		
		if (!flip(surv)) {
			do death;
		}
	}
	
	reflex activity_choice {
		
		if (fed) {
			hours_since_feed <- hours_since_feed + 1;
		}
		
		if (num_eggs > 0) {
			hours_since_ovi <- hours_since_ovi + 1;
		}
		
		
//		If the female has mated and has blood fed, she can start producing eggs
		if (!virgin and !gravid and fed) {
			gravid <- true;
		}
		
//		If she hasn't previously fed, but is old enough to, she will get hungry
		if (!has_bitten and cycle >= age_hours + hours_feed) {
			hungry <- true;
		}
//		If she's already gone through a gonotrophic cycle, isn't currently gravid, and enough time as past, she will get hungry
		if (num_eggs > 0 and !gravid and hours_since_ovi >= hours_refeed) {
			hungry <- true;
		}
		
//		Incubation period to produce eggs for oviposition
		if (gravid and egg_timer < hours_ovi) {
			egg_timer <- egg_timer + 1;
		}
		
//		Decide activity based on time of day
		if (morning or afternoon) {
			if (!mated) {
				do Mate;
			}
			
			else if (gravid and egg_timer >= hours_ovi) {
				do Lay_eggs;
			}

			else if (hungry) {
				do Feed;
			}
		}
	}
	
	
	action Lay_eggs {
		create (egg) number: ovi_eggs {
			if (myself.fsRIDL) {	
				self.fsRIDL_carrier <- true;
			}
			if (myself.fsRIDL_carrier) {
				if(flip(0.5)) {
					self.fsRIDL_carrier <- true;	
				}
			}
			myself.num_eggs <- myself.num_eggs + 1;
		}
		egg_timer <- 0;
		gravid <- false;
		fed <- false;
		hours_since_ovi <- 0;
		num_GC <- num_GC + 1;
		GC_plus_len <- GC_plus_len + 1;
	}
	
	action Feed {
		if (!flip(biting_success)) {
			num_bites <- num_bites + 1;
			hungry <- false;
			fed <- true;
			hours_since_feed <- 0;
			
			if (nb_days >= release_delay) {
				bite_count <- bite_count + 1;
				if (has_bitten) {
					bite_plus_count <- bite_plus_count + 1;
				}
			}
			if (!has_bitten) {
				has_bitten <- true;
			}
		}
	}
	
	action death {
		float daily_bite <- 0.0;
		float GC_bite <- num_bites / num_GC;
		
		if (num_GC > 1) {
			GC_plus_len <- GC_plus_len - (num_GC - 1);
		}
		
		if (age > 0) {
			daily_bite <- num_bites / age;
		}
		if (nb_days >= release_delay and is_batch) {
			save [	self.name,
					cycle,
					exp_TMT,
					exp_SIT,
					exp_fsRIDL,
					exp_TMT_reduced,
					age,
					num_mates,
					num_bites,
					num_GC,
					self.death_cause,
					daily_bite,
					GC_bite
				]
			to: "TMT_female_output.csv" format: csv rewrite: false;
		}
		do die;
	}
	
	action Mate{
//		The female will randomly chose one male from the total pool of WT and Tx males
//		She will be less likely to choose a transgenic male depending on their competitiveness
		int male_choice <- rnd(length(male_mosquito) + round((length(transgenic_mosquito) * tx_competitiveness)));
		if (male_choice <= length(male_mosquito) and length(male_mosquito) > 0) {
			male_mosquito mate <- one_of(male_mosquito);
			ask mate {
				if (!flip(mating_prob)) {
					myself.mated <- true;
					myself.virgin <- false;
					
//					Re-mating with WT male displaces fsRIDL sperm
					if (myself.fsRIDL) {
						myself.fsRIDL <- false;
					}
//					If the male is a carrier for fsRIDL, half of their offspring will be too
					if (self.fsRIDL_carrier) {
						myself.fsRIDL_carrier <- true;					
					}
//					If this is my first mate, begin refractory period
					if (myself.mated_timer = 0) {
						myself.mated_timer <- 1;
					}
					
					myself.num_mates <- myself.num_mates + 1;
					
					if (myself.num_mates = 1) {
						mate1 <- mate1 + 1;
					}
					if (myself.num_mates = 2) {
						mate1 <- mate1 - 1;
						mate2 <- mate2 + 1;
					}
					if (myself.num_mates = 3) {
						mate2 <- mate2 - 1;
						mate3 <- mate3 + 1;
					}
					if (myself.num_mates = 4) {
						mate3 <- mate3 - 1;
						mateplus <- mateplus + 1;
					}
				}
			}
		}
		else if (length(transgenic_mosquito) > 0) {
			transgenic_mosquito mate <- one_of(transgenic_mosquito);
			ask mate {
				if (!flip(mating_prob)) {
					myself.mated <- true;
					myself.virgin <- false;
				}
				if (myself.mated and fsRIDL) {
					myself.fsRIDL <- true;
//					Re-mating with a fsRIDL male displaces fsRIDL carrier sperm
					if (myself.fsRIDL_carrier) {
						myself.fsRIDL_carrier <- false;
					}
				}
				if (myself.mated and TMT) {
					if (flip(TMT_lethality)) {
						myself.toxin <- true;	
					}
				}
				if (myself.mated and SIT) {
					myself.virgin <- true;
				}
			}
			if (toxin) {
				self.death_cause <- "TMT";
				do death;
			}
		}
	}
}

species male_mosquito parent: mosquito {
	
	
	reflex age when: new_day{
		age <- age + 1;
		
		float surv <- (-male_delta * (age)) + (1 + male_delta);
		
		if (!flip(surv)) {
			do die;
		}
	}
}

species transgenic_mosquito parent: male_mosquito {
	
	bool TMT <- false;

	init {
		if (exp_fsRIDL = true) {
			self.fsRIDL <- true;
		}
		else if (exp_TMT = true) {
			self.TMT <- true;
		}
		else if (exp_TMT_reduced = true) {
			self.TMT <- true;
		}
		else if (exp_SIT = true) {
			self.SIT <- true;
		}
	}
}

experiment TMT_batch type: batch repeat: 10 keep_seed: false until: nb_days >= (release_delay + 120) {
	parameter "TMT" var:exp_TMT init: true;
	parameter "Batch" var:is_batch init: true;
	
	output {
		monitor Nb_days value: nb_days refresh: every(day #cycles);
		monitor Days_since_release value: days_since_release refresh: every(day #cycles);
	}
}

experiment TMT_reduced_batch type: batch repeat: 10 keep_seed: false until: nb_days >= (release_delay + 120) {
	parameter "TMT_reduced" var:exp_TMT_reduced init: true;
	parameter "Batch" var:is_batch init: true;
	parameter "Lethality" var:TMT_lethality init: TMT_reduced_lethality;
	
	output {
		monitor Nb_days value: nb_days refresh: every(day #cycles);
		monitor Days_since_release value: days_since_release refresh: every(day #cycles);
	}
}

experiment fsRIDL_batch type: batch repeat: 10 keep_seed: false until: nb_days >= (release_delay + 120) {
	parameter "fsRIDL" var:exp_fsRIDL init: true;
	parameter "Batch" var:is_batch init: true;
	
	output {
		monitor Nb_days value: nb_days refresh: every(day #cycles);
		monitor Days_since_release value: days_since_release refresh: every(day #cycles);
	}
}


experiment SIT_batch type: batch repeat: 10 keep_seed: false until: nb_days >= (release_delay + 120) {
	parameter "SIT" var:exp_SIT init: true;
	parameter "Batch" var:is_batch init: true;
	
	output {
		monitor Nb_days value: nb_days refresh: every(day #cycles);
		monitor Days_since_release value: days_since_release refresh: every(day #cycles);
	}
}
