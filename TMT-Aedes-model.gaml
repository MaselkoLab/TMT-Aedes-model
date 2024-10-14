/**
* Name: TMT-Aedes-model
* Author: Samuel Beach 
* Version: 3.1.0
* Date: 2024-10-04
*/


model TMT

global {
	
//	Experimental variables

	float tx_competitiveness <- 0.031;	
	int release_period <- 3;
	int release_delay <- 150;				
	float TMT_lethality <- 0.75;
	
	float release_ratio <- 12.0;			// Set to 3.0, 6.0, or 12.0
	float polyandry <- 0.2;					// Set to 0.0, 0.2, or 0.55
	float DDM <- 0.0075;					// Set to 0.0025, 0.0075, or 0.025

	float DDM_x;
	float DDM_a;
	float DDM_m;
	float DDM_k;
	
	float polyandry_m;
	float polyandry_b;
	
//	Global time variables
	int day <- 24;
	int week <- day * 7;
	int year <- day * 365;
	
//	Initialisation parameters
	int nb_init <- 10000;
	bool is_batch <- false;
	
//	Immature stage-specific parameters
	float imm_surv <- 0.9975;
	int nb_larvae function: length(larvae) > 0? length(larvae) + length(larvae_L3): 0;
	
	float DDM_surv function: nb_larvae > 0? (DDM_x - DDM_a) / (1 + exp(DDM_k * (nb_larvae - DDM_m))) + DDM_a : DDM_x;
	
//	Mosquito-specific parameters
	float biting_success <- 0.8;
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
	bool exp_TMT <- false;
	bool exp_SIT <- false;
	bool exp_WT <- false;
	bool exp_TMT_SIT <- false;
	string condition <- "WT";
	
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
	int nb_males_init;
	int nb_moz_init;
	int nb_moz_release;
	int max_moz_count;
	int daily_moz;
	int total_moz;
	int mean_moz <- 1;
	int days_since_release function: nb_days >= release_delay?
										nb_days - release_delay:
										0;
	
	reflex max_count when: (length(male_mosquito) + length(female_mosquito)) > max_moz_count {
		max_moz_count <- (length(male_mosquito) + length(female_mosquito));
	}
	
	reflex release when: !exp_WT and
						 new_day and
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
				self.cycle,
				exp_WT,
				exp_TMT,
				exp_fsRIDL,
				exp_SIT,
				exp_TMT_SIT,
				DDM,
				polyandry,
				release_ratio,
				TMT_lethality,
				p50,
				self.bite_count,
				self.bite_plus_count
		]
		to: "TMT_PR50.csv" format: csv rewrite: false;
	}
	
//	Time to 95% female population reduction
	reflex when: nb_days >= release_delay and p95 = 0 and length(female_mosquito) <= nb_females_init * 0.05 {
		p95 <- cycle - (release_delay*24);
		save [	self.name,
				self.cycle,
				exp_WT,
				exp_TMT,
				exp_fsRIDL,
				exp_SIT,
				exp_TMT_SIT,
				DDM,
				polyandry,
				release_ratio,
				TMT_lethality,
				p95,
				self.bite_count,
				self.bite_plus_count
		]
		to: "TMT_PR95.csv" format: csv rewrite: false;
	}
	
	reflex step_save when: is_batch and nb_days >= release_delay  {
		save [	self.name,
				self.cycle,
				length(female_mosquito),
				length(male_mosquito),
				self.nb_females_init,
				self.nb_males_init,
				self.nb_moz_init,
				GC_plus_len,
				GC_plus_init,
				self.bite_count,
				self.bite_plus_count,
				self.sum_Tx,
				exp_WT,
				exp_TMT,
				exp_fsRIDL,
				exp_SIT,
				exp_TMT_SIT,
				DDM,
				polyandry,
				release_period,
				release_ratio,
				TMT_lethality,
				tx_competitiveness
		]
		to: "TMT_step_output.csv" format: csv rewrite: false;
	}
	
	reflex init_female_count when: cycle = (release_delay * 24) {
		nb_females_init <- length(female_mosquito);
		nb_males_init <- length(male_mosquito);
		nb_moz_init <- (length(female_mosquito) + length(male_mosquito));
		GC_plus_init <- GC_plus_len;
	}

	reflex end when: (nb_days > release_delay) and length(female_mosquito) = 0 and !is_batch {
		do pause;
	}
	
	init {
		
//		Initialise the values determining the strength of density dependent mortality
// 		y = (x-a)/(1+e^(k(n-m)))+a
//		y = adjusted hourly survival rate
//		x = Upper bound survival rate
//		a = Lower bound survival rate
//		n = Number of larvae
//		m = Midpoint of the sigmoidal curve (equilibrium point)
//		k = Steepness of the curve

		DDM_x <- 0.97 + DDM;
		DDM_a <- 0.97 - DDM;
		DDM_m <- 300.0;
		DDM_k <- 0.0001 / DDM;
		
		
		if (polyandry = 0.0) {			// remating =  0%
			polyandry_m <- -0.0;	
			polyandry_b <- 0.0;	 
		}
		else if (polyandry = 0.2) {		// remating ≈ 20%
			polyandry_m <- -0.16;
			polyandry_b <- 0.16;
		}
		else if (polyandry = 0.55) {		// remating ≈ 55%
			polyandry_m <- -0.0475;
			polyandry_b <- 0.14;
		}
		
		
		if (condition = "WT") {
			exp_WT <- true;
		}
		else if (condition = "fsRIDL") {
			exp_fsRIDL <- true;
		}
		else if (condition = "SIT") {
			exp_SIT <- true;
		}
		else if (condition = "TMT") {
			exp_TMT <- true;
		}
		else if (condition = "TMT_SIT") {
			exp_TMT_SIT <- true;
		}
		
		create egg number: nb_init {
			age <- round(gauss(96,12));	
		}
		create female_mosquito number: 150 {
			
		}
		create male_mosquito number: 150 {
			
		}
	}
}

species egg  {
	int age;
	int hours_hatch;
	bool fsRIDL_carrier;
	bool male;
	bool female;
	float surv;
	string sperm;
	
	init {
		age <- cycle;
		hours_hatch <- round(gauss(96,12));
		if (flip(0.5)) {
			female <- true;
		}
		else {
			male <- true;
		}
		surv <- DDM_surv;
		
	}
	
	reflex age {		
		if (sperm = "SIT") {
			do die;
		}
		if (DDM_surv < surv) {
			surv <- DDM_surv;
		}
		if (!flip(surv)) {
			do die;
		}
		if (cycle >= age + hours_hatch) {
			create (larvae) {
				self.fsRIDL_carrier <- myself.fsRIDL_carrier;
				self.male <- myself.male;
				self.female <- myself.female;
			}
			do die;
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
		surv <- DDM_surv;
	}

	reflex age {
		if (DDM_surv < surv) {
			surv <- DDM_surv;
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
			if (female) {
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
	int age_hours;
	bool SIT <- false;
	bool fsRIDL <- false;
	bool fsRIDL_carrier <- false;
	bool receptive <- false;
	
	init {
		age_hours <- cycle;
	}
	
//	Mosquitoes will beceome receptive to mating 24 - 48 hours after eclosion
	reflex receptive when: !receptive and (cycle - age_hours >= 24) {		
		if (flip(1/24)) {
			receptive <- true;
		}
	}
}

species female_mosquito parent: mosquito {
	bool gravid <- false;
	bool mated <- false;
	bool fed <- false;
	bool hungry <- false;
	float refeeding <- 0.0003;
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
	float remating;
	
	list<string> spermatheca;
	
	int ovi_eggs;
	bool has_bitten;
	
	init {
		hours_feed <- round(gauss(36,6));
		hours_ovi <- round(gauss(72,12));
		hours_refeed <- round(gauss(12,4));
		ovi_eggs <- round(gauss(60,10));
	}
	
	reflex polyandry when: mated and (mated_timer > 0 and mated_timer < 25) {
		remating <- (((polyandry_m)*ln(mated_timer)) + polyandry_b);
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
//		if (!virgin and !gravid and fed) {
		if (!gravid and fed) {
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
		
//		Activity selection
		if (morning or afternoon) {
			if (!mated and receptive) {
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
			
			sperm <- one_of(myself.spermatheca);
			
			if (sperm = "fsRIDL") {
				self.fsRIDL_carrier <- true;
			}
			else if (sperm = "fsRIDL_carrier") {
				if (flip(0.5)) {
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
		if (flip(biting_success)) {	
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
					exp_WT,
					exp_TMT,
					exp_SIT,
					exp_fsRIDL,
					exp_TMT_SIT,
					DDM,
					polyandry,
					release_ratio,
					TMT_lethality,
					age,
					num_mates,
					num_bites,
					num_eggs,
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
//		The female will randomly chose one male from the total pool of receptive WT and Tx males
//		She will be less likely to choose a transgenic male depending on their competitiveness
		int male_choice <- rnd(length(male_mosquito where (each.receptive)) + round((length(transgenic_mosquito where (each.receptive)) * tx_competitiveness)));
		if (male_choice <= length(male_mosquito where (each.receptive)) and length(male_mosquito where (each.receptive)) > 0) {
			male_mosquito mate <- one_of(male_mosquito where receptive);
			ask mate {
				if (flip(mating_prob)) {
					myself.mated <- true;
					
					if (self.fsRIDL_carrier) {
						add "fsRIDL_carrier" to: myself.spermatheca;
					}
					else {
						add "WT" to: myself.spermatheca;
					}

//					If this is my first mate, begin refractory period
					if (myself.mated_timer = 0) {
						myself.mated_timer <- 1;
					}
					
					myself.num_mates <- myself.num_mates + 1;
					
					if (!is_batch) {
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
		}
		else if (length(transgenic_mosquito) > 0) {
			transgenic_mosquito mate <- one_of(transgenic_mosquito where receptive);
			ask mate {
				if (flip(mating_prob)) {
					myself.mated <- true;
//					If this is my first mate, begin refractory period
					if (myself.mated_timer = 0) {
						myself.mated_timer <- 1;
					}
					
					if (self.fsRIDL) {
						add "fsRIDL" to: myself.spermatheca;
					}
					if (self.SIT) {
						add "SIT" to: myself.spermatheca;
					}
					if (self.TMT) {
						if (flip(TMT_lethality)) {
							myself.toxin <- true;	
						}
					}
					
					myself.num_mates <- myself.num_mates + 1;
					
					if (!is_batch) {
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
		else if (exp_SIT = true) {
			self.SIT <- true;
		}
		else if (exp_TMT_SIT) {
			self.SIT <-true;
			self.TMT <- true;
		}
	}
}


experiment WT type: gui {
	
	parameter "WT" var:exp_WT init: true;
	
	output {
		monitor Nb_days value: nb_days refresh: every(day #cycles);
		monitor Nb_mosquitoes value: (length(female_mosquito)+length(male_mosquito)) refresh: every(day #cycles);
		monitor Nb_agents value: (length(egg)+length(larvae)+length(larvae_L3)+length(pupae)+length(female_mosquito)+length(male_mosquito)) refresh: every(day #cycles);
		monitor Nb_L1L4 value: nb_larvae refresh: every(day #cycles);
		monitor DDM_surv value: DDM_surv refresh: every(day #cycles);
		monitor Polyandry value: polyandry;
		monitor DDM value: DDM;
		
		
		display nbMoz type: java2D refresh: every(day #cycles) {
			chart "Wild mosquito population numbers" type: series {
				data "Females" value: length(female_mosquito) color: #red;
				data "Males" value: length(male_mosquito) color: #blue;
				data "Total" value: (length(male_mosquito) + length(female_mosquito)) color: #black;
			}
		}
		
		display DDM type: java2D refresh: every(day #cycles) {
			chart "Larval survival rate" type: series y_range: [0.967,1.0] {
				data "Hourly survival rate" value: (DDM_surv) ;
			}
		}

		
		display MateCount type: java2D refresh: every(day #cycles) {
			chart "How many mates has each female had??" type: pie {
				data "1" value: mate1 color: #red;
				data "2" value: mate2 color: #blue;
				data "3" value: mate3 color: #green;
				data "4+" value: mateplus color: #yellow;
			}
		}
		
		display Proportion type: java2D refresh: every(day #cycles) {
			chart "Current mosquito population breakdown" type: pie {
				data "Wild Females" value: length(female_mosquito) color: #red;
				data "Wild Males" value: length(male_mosquito) color: #blue;
				data "Transgenic Males" value: length(transgenic_mosquito) color: #green;
			}
		}  
	}
}

experiment TMT type: gui {
	
	parameter "TMT" var:condition init: "TMT";
	
	output {
		monitor Nb_days value: nb_days refresh: every(day #cycles);
		monitor Nb_mosquitoes value: (length(female_mosquito)+length(male_mosquito)) refresh: every(day #cycles);
		monitor Nb_agents value: (length(egg)+length(larvae)+length(larvae_L3)+length(pupae)+length(female_mosquito)+length(male_mosquito)) refresh: every(day #cycles);
		monitor Nb_L1L4 value: nb_larvae refresh: every(day #cycles);
		monitor DDM_surv value: DDM_surv refresh: every(day #cycles);
		monitor Polyandry value: polyandry;
		monitor DDM value: DDM;
		monitor Release_ratio value: release_ratio;
		monitor Release_period value: release_period;
		monitor Total_nb_transgenics value: sum_Tx refresh: every(day #cycles);
		monitor Bite_count value: bite_count refresh: every(day #cycles);
		monitor Bite_plus_count value: bite_plus_count refresh: every(day #cycles);
		monitor p50 value: p50 refresh: every(day #cycles);
		monitor p95 value: p95 refresh: every(day #cycles);
		
		display nbMoz type: java2D refresh: every(day #cycles) {
			chart "Wild mosquito population numbers" type: series {
				data "Females" value: length(female_mosquito) color: #red;
				data "Males" value: length(male_mosquito) color: #blue;
				data "Total" value: (length(male_mosquito) + length(female_mosquito)) color: #black;
			}
		}
		
		display DDM type: java2D refresh: every(day #cycles) {
			chart "Larval survival rate" type: series y_range: [0.967,1.0] {
				data "Hourly survival rate" value: (DDM_surv) ;
			}
		}

		
		display MateCount type: java2D refresh: every(day #cycles) {
			chart "How many mates has each female had??" type: pie {
				data "1" value: mate1 color: #red;
				data "2" value: mate2 color: #blue;
				data "3" value: mate3 color: #green;
				data "4+" value: mateplus color: #yellow;
			}
		}
		
		display Proportion type: java2D refresh: every(day #cycles) {
			chart "Current mosquito population breakdown" type: pie {
				data "Wild Females" value: length(female_mosquito) color: #red;
				data "Wild Males" value: length(male_mosquito) color: #blue;
				data "Transgenic Males" value: length(transgenic_mosquito) color: #green;
			}
		}  
	}
}

experiment TMT_SIT type: gui {
	
	parameter "TMT_SIT" var:condition init: "TMT_SIT";
	
	output {
		monitor Nb_days value: nb_days refresh: every(day #cycles);
		monitor Nb_mosquitoes value: (length(female_mosquito)+length(male_mosquito)) refresh: every(day #cycles);
		monitor Nb_agents value: (length(egg)+length(larvae)+length(larvae_L3)+length(pupae)+length(female_mosquito)+length(male_mosquito)) refresh: every(day #cycles);
		monitor Nb_L1L4 value: nb_larvae refresh: every(day #cycles);
		monitor DDM_surv value: DDM_surv refresh: every(day #cycles);
		monitor Polyandry value: polyandry;
		monitor DDM value: DDM;
		monitor Release_ratio value: release_ratio;
		monitor Release_period value: release_period;
		monitor Total_nb_transgenics value: sum_Tx refresh: every(day #cycles);
		monitor Bite_count value: bite_count refresh: every(day #cycles);
		monitor Bite_plus_count value: bite_plus_count refresh: every(day #cycles);
		monitor p50 value: p50 refresh: every(day #cycles);
		monitor p95 value: p95 refresh: every(day #cycles);
		
		display nbMoz type: java2D refresh: every(day #cycles) {
			chart "Wild mosquito population numbers" type: series {
				data "Females" value: length(female_mosquito) color: #red;
				data "Males" value: length(male_mosquito) color: #blue;
				data "Total" value: (length(male_mosquito) + length(female_mosquito)) color: #black;
			}
		}
		
		display DDM type: java2D refresh: every(day #cycles) {
			chart "Larval survival rate" type: series y_range: [0.967,1.0] {
				data "Hourly survival rate" value: (DDM_surv) ;
			}
		}

		
		display MateCount type: java2D refresh: every(day #cycles) {
			chart "How many mates has each female had??" type: pie {
				data "1" value: mate1 color: #red;
				data "2" value: mate2 color: #blue;
				data "3" value: mate3 color: #green;
				data "4+" value: mateplus color: #yellow;
			}
		}
		
		display Proportion type: java2D refresh: every(day #cycles) {
			chart "Current mosquito population breakdown" type: pie {
				data "Wild Females" value: length(female_mosquito) color: #red;
				data "Wild Males" value: length(male_mosquito) color: #blue;
				data "Transgenic Males" value: length(transgenic_mosquito) color: #green;
			}
		}  
	}
}

experiment fsRIDL type: gui {
	
	parameter "fsRIDL" var:condition init: "fsRIDL";
	
	output {
		monitor Nb_days value: nb_days refresh: every(day #cycles);
		monitor Nb_mosquitoes value: (length(female_mosquito)+length(male_mosquito)) refresh: every(day #cycles);
		monitor Nb_agents value: (length(egg)+length(larvae)+length(larvae_L3)+length(pupae)+length(female_mosquito)+length(male_mosquito)) refresh: every(day #cycles);
		monitor Nb_L1L4 value: nb_larvae refresh: every(day #cycles);
		monitor DDM_surv value: DDM_surv refresh: every(day #cycles);
		monitor Polyandry value: polyandry;
		monitor DDM value: DDM;
		monitor Release_ratio value: release_ratio;
		monitor Release_period value: release_period;
		monitor Total_nb_transgenics value: sum_Tx refresh: every(day #cycles);
		monitor Bite_count value: bite_count refresh: every(day #cycles);
		monitor Bite_plus_count value: bite_plus_count refresh: every(day #cycles);
		monitor p50 value: p50 refresh: every(day #cycles);
		monitor p95 value: p95 refresh: every(day #cycles);
		
		display nbMoz type: java2D refresh: every(day #cycles) {
			chart "Wild mosquito population numbers" type: series {
				data "Females" value: length(female_mosquito) color: #red;
				data "Males" value: length(male_mosquito) color: #blue;
				data "Total" value: (length(male_mosquito) + length(female_mosquito)) color: #black;
			}
		}
		
		display DDM type: java2D refresh: every(day #cycles) {
			chart "Larval survival rate" type: series y_range: [0.967,1.0] {
				data "Hourly survival rate" value: (DDM_surv) ;
			}
		}

		
		display MateCount type: java2D refresh: every(day #cycles) {
			chart "How many mates has each female had??" type: pie {
				data "1" value: mate1 color: #red;
				data "2" value: mate2 color: #blue;
				data "3" value: mate3 color: #green;
				data "4+" value: mateplus color: #yellow;
			}
		}
		
		display Proportion type: java2D refresh: every(day #cycles) {
			chart "Current mosquito population breakdown" type: pie {
				data "Wild Females" value: length(female_mosquito) color: #red;
				data "Wild Males" value: length(male_mosquito) color: #blue;
				data "Transgenic Males" value: length(transgenic_mosquito) color: #green;
			}
		} 
	}
}

experiment SIT type: gui {
	
	parameter "SIT" var:condition init: "SIT";
	
	output {
		monitor Nb_days value: nb_days refresh: every(day #cycles);
		monitor Nb_mosquitoes value: (length(female_mosquito)+length(male_mosquito)) refresh: every(day #cycles);
		monitor Nb_agents value: (length(egg)+length(larvae)+length(larvae_L3)+length(pupae)+length(female_mosquito)+length(male_mosquito)) refresh: every(day #cycles);
		monitor Nb_L1L4 value: nb_larvae refresh: every(day #cycles);
		monitor DDM_surv value: DDM_surv refresh: every(day #cycles);
		monitor Polyandry value: polyandry;
		monitor DDM value: DDM;
		monitor Release_ratio value: release_ratio;
		monitor Release_period value: release_period;
		monitor Total_nb_transgenics value: sum_Tx refresh: every(day #cycles);
		monitor Bite_count value: bite_count refresh: every(day #cycles);
		monitor Bite_plus_count value: bite_plus_count refresh: every(day #cycles);
		monitor p50 value: p50 refresh: every(day #cycles);
		monitor p95 value: p95 refresh: every(day #cycles);
		
		display nbMoz type: java2D refresh: every(day #cycles) {
			chart "Wild mosquito population numbers" type: series {
				data "Females" value: length(female_mosquito) color: #red;
				data "Males" value: length(male_mosquito) color: #blue;
				data "Total" value: (length(male_mosquito) + length(female_mosquito)) color: #black;
			}
		}
		
		display DDM type: java2D refresh: every(day #cycles) {
			chart "Larval survival rate" type: series y_range: [0.967,1.0] {
				data "Hourly survival rate" value: (DDM_surv) ;
			}
		}

		
		display MateCount type: java2D refresh: every(day #cycles) {
			chart "How many mates has each female had??" type: pie {
				data "1" value: mate1 color: #red;
				data "2" value: mate2 color: #blue;
				data "3" value: mate3 color: #green;
				data "4+" value: mateplus color: #yellow;
			}
		}
		
		display Proportion type: java2D refresh: every(day #cycles) {
			chart "Current mosquito population breakdown" type: pie {
				data "Wild Females" value: length(female_mosquito) color: #red;
				data "Wild Males" value: length(male_mosquito) color: #blue;
				data "Transgenic Males" value: length(transgenic_mosquito) color: #green;
			}
		} 
	}
}

experiment bulk_batch type: batch repeat: 10 keep_seed: false until: nb_days >= (release_delay + 300) {
	parameter "Batch" var:is_batch init: true;

	parameter "DDM" var:DDM among: [0.0025, 0.0075, 0.025];
	parameter "Polyandry" var:polyandry among: [0.0, 0.2, 0.55];
	parameter "Release ratio" var:release_ratio among: [3.0, 6.0, 12.0];
	parameter "Biocontrol method" var:condition among: ["WT", "fsRIDL", "SIT", "TMT"];
	
	output {
		monitor Nb_days value: nb_days refresh: every(day #cycles);
	}
}

experiment lethality_batch type: batch repeat: 10 keep_seed: false until: nb_days >= (release_delay + 300) {
	parameter "Batch" var:is_batch init: true;

	parameter "DDM" var:DDM init: 0.0075;
	parameter "Polyandry" var:polyandry init: 0.2;
	parameter "Release ratio" var:release_ratio init: 12.0;
	
	parameter "TMT sterility" var:condition among: ["TMT", "TMT_SIT"];
	parameter "TMT lethality" var:TMT_lethality min: 0.1 max: 1.0 step: 0.1;
	
	output {
		monitor Nb_days value: nb_days refresh: every(day #cycles);
	}
}
