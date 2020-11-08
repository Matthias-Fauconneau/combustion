//let [(_, νf), (_, νr)] = &equation;
//let [Σνf, Σνr] = [νf.iter().sum::<u8>() as f64, νr.iter().sum::<u8>() as f64];
/*fn dT_Cp(&self, T: f64) -> f64 {
	use itertools::Itertools;
	let a = &self.coefficients[self.temperature_ranges.iter().tuple_windows().position(|(&min, &max)| min <= T && T <= max).unwrap_or_else(|| panic!("{:?}", T))];
	ideal_gas_constant * (a[1]+2.*a[2]*T+3.*a[3]*T*T+4.*a[4]*T*T*T)
}
fn dT_b(&self, T: f64) -> f64 { // dT(S/R - H/RT)
	use itertools::Itertools;
	let a = &self.coefficients[self.temperature_ranges.iter().tuple_windows().position(|(&min, &max)| min <= T && T <= max).unwrap_or_else(|| panic!("{:?}", T))];
	(a[0]-1.)/T + a[1]/2. + a[2]/12.*T + a[3]/36.*T*T + a[4]/80.*T*T*T + a[5]/(T*T)
}*/
/*impl<const S: usize> State<S> {
	pub fn step(&mut self, System{thermodynamics: species, reactions, amount, volume: V, molar_masses: W}: &System<S>) {
		use iter::array::{from_iter, map, generate};
		let scale = |s, v| from_iter(scale(s, v.iter().copied()));
		macro_rules! map { ($($args:ident),*| $expr:expr) => { izip!($($args),*).map(|($($args),*)| $expr) } }

		let rcpn = 1. / amount;
		let T = self.temperature;
		let B = map(species, |s| s.b(T));
		let dTB = map(species, |s| s.dT_b(T));
		let rcpV = 1. / V;
		// C = n / V
		let C = scale(rcpV, self.amounts);
		let a = S-1;
		let Ca = C[a];
		let mut ω = vec!(0.; S-1); // Skips most abundant specie (last index) (will be deduced from conservation)
		let mut dTω = vec!(0.; S-1);
		let mut dVω = vec!(0.; S-1);
		let mut dnω = vec!(vec!(0.; S-1); S-1); //[k][j]
		for Reaction{equation, rate_constant: rate_constant@RateConstant{temperature_exponent: β, activation_energy: Ea}, model, specie_net_coefficients: ν, Σνf, Σνr, PRν} in reactions.iter() {
			let equilibrium_constant = PRν * f64::exp(dot(ν, B));
			let kf = arrhenius(rate_constant, T);
			let kr = kf / equilibrium_constant;
			// ΠC^ν
			let [ΠCνf, ΠCνr] = from_iter(equation.iter().map(|(species, ν)| species.iter().zip(ν.iter()).map(|(&specie, &ν)| C[specie].powi(ν as i32)).product::<f64>() ));
			let Rf = kf * ΠCνf;
			let Rr = kr * ΠCνr;
			let νfRfνrRr = vec!(0.; S);
			let [forward, reverse] = equation;
			for (specie, ν) in forward { νfRfνrRr[specie] += ν * Rf; }
			for (specie, ν) in reverse { νfRfνrRr[specie] -= ν * Rr; }
			let c = model.efficiency(T, &C, kf);
			let R = Rf - Rr;
			// dT(R) = (β+Ea/(T.Rideal))/T.R + Rr. Σ ν.dT(B) - νfRfνrRr[a]/T
			let dTR = (β+Ea/(T*ideal_gas_constant))/T*R + Rr*dot(ν,dTB) - νfRfνrRr[a] / T;
			// dV(R) = 1/V . ( (kf.Sf-kr.Sr) - (Σνf.Rf - Σνr.Rr) )
			let dVR = rcpV * ( νfRfνrRr[a] - (Σνf*Rf - Σνr*Rr));
			// dn(R) = 1/n . ( kf.(Sf-Sfa) - kr.(Sr-Sra) )
			let dnR = map(νfRfνrRr, |νfRfνrRrj| rcpn * (νfRfνrRrj - νfRfνrRr[a]));
			let (dTc, dVc, dnc) = match model {
				Model::Elementary => (0., 0., vec!(0.; S-1)),
				Model::ThirdBody{efficiencies}|Model::Falloff{efficiencies} => {
					let has = map(efficiencies, |e| if e != 0. { 1. } else { 0. });
					(
						// dT(c) = has(a) . (-C/T)
						has[a] * -C/T,
						// dV(c) = 1/V . ( has(a). C - Σ C )
						rcpV * (has[a] * C  - dot(has, C)),
						// dn(c) = 1/V . ( has(n) - has(a) )
						has[..a].map(|has_n| rcpV * (has_n - has[a]))
					)
				}
			};
			let cR = c * R;
			let RdTccdTR = R * dTc + c * dTR;
			let RdVccdVR = R * dVc + c * dVR;
			let RdnccdnR = from_iter(map!(dnc,dnR| R*dnc + c*dnR));
			for (specie, ν) in ν.iter().enumerate() {
				//ω̇̇̇̇̇̇̇̇̇̇ = Σ ν c R
				ω[specie] += ν * cR;
				// dT(ω̇̇̇̇̇̇̇̇̇̇) = Σ ν.(R.dT(c)+c.dT(R))
				dTω[specie] += ν * RdTccdTR;
				// dV(ω) = Σ ν.(R.dV(c)+c.dV(R))
				dVω[specie] += ν * RdVccdVR;
				// dn(ω) = Σ ν.(R.dn(c)+c.dn(R))
				for dnω in dnω[specie] { dnω += ν * RdnccdnR; }
			}
		}
		use nalgebra::{Matrix, MatrixMN};
		let mut J = unsafe{MatrixMN::<f64, {2+S-1}, {2+S-1}>::new_uninitialized()}; // fixme
		let Cp = map(species, |s| s.specific_heat_capacity(T));
		// 1 / (Σ C.Cp)
		let rcp_ΣCCp = 1./dot(C, Cp);
		let H = species.map(|s| s.specific_enthalpy(T));
		let Wa = W[a];
		let HaWa = H[a]/Wa;
		// Ha/Wa*W - H
		let HaWaWH = from_iter(map!(W,H| HaWa*W - H));
		// dtT = - 1 / (Σ C.Cp) . Σ H.ω̇
		let dtT = - rcp_ΣCCp * dot(H, ω);
		let CpaWa = Cp[a]/Wa;
		// dT(dtT) = 1 / (Σ C.Cp) . [ dtT . Σ C.(Cpa/T - dT(Cp)) + Σ_ ( (Ha/Wa*W - H).dT(ω) + (Cpa/Wa.W - Cp).ω ) ]
		let dTdtT = rcp_ΣCCp * ( dtT * dot(C, species.iter().map(|s| Cp[a]/T - s.dt_Cp())) + dot(HaWaWH, dTω) + dotia(map!(W,Cp, CpaWa*W - Cp), ω));
		J[(0,0)] = dTdtT;
		// dV(dtT) = 1 / (Σ C.Cp) . [ Σ_ (Ha/Wa*W - H).dV(ω) + dtT/V . Σ_ C.(Cp-Cpa) ]
		let dVdtT = rcp_ΣCCp * ( dot(HaWaWH, dVω) + rcpV * dtT * dotai(C, Cp.iter().map(|Cp| Cp - Cp[a])));
		J[(0,1)] = dVdtT;
		// dn(dtT) = 1 / (Σ C.Cp) . [ Σ_ (Ha/Wa*W - H).dn(ω) + dtT/V . (Cpa-Cp) ]
		let dndtT = from_iter(map!(dnω, Cp| rcp_ΣCCp * ( dot(HaWaWH, dnω) + rcpV * dtT * (Cp[a]-Cp) )));
		J.row_part_mut(0,2,S+1).copy_from_slice(dndtT);

		// dT(dtV) = V/C . Σ_ (1-W/Wa).(dT(ω)+ω/T) + V/T.(dT(dtT) - dtT/T)
		let rcpC = V*rcpn;
		let WWa = map(W, |W| 1.-W/Wa);
		let VT = V/T;
		let dTdtV = V*rcpC * map!(WWa,dTω,ω| WWa*(dTω+ω/T)).sum() + VT * (dTdtT - dtT/T);
		J[(1,0)] = dTdtV;
		// dV(dtn) = VdV(ω)+ω
		let dVdtn = from_iter(map!(dVω,ω| V*dVω+ω));
		// dV(dtV) = 1/C . Σ_ (1-W/Wa).dV(dtn) + 1/T.(V.dV(dtT)+dtT)
		let dVdtV = rcpC * dot(WWa, dVdtn) + 1./T*(V*dVdtT+dtT);
		J[(1,1)] = dVdtV;
		// dn(dtn) = Vdn(ω)
		let dndtn = generate(S-1, |j| generate(S-1, |k| V*dnω[k][j])); // Transpose [k][j] -> [j][k]
		// dn(dtV) = 1/C . Σ_ (1-W/Wa).dn(dtn) + V/T.dn(dtT))
		let dndtV = map!(dndtn, dndtT| rcpC * dot(WWa, dndtn)+ VT*dndtT);
		J.row_part_mut(1,2,S+1) = Matrix::from_iterator(dndtV);

		// dT(dtn) = VdT(ω)
		let dTdtn = scale(V, dTω);
		J.column_part_mut(0,2,S+1) = Matrix::from_iterator(dTdtn);
		// dV(dtn)
		J.column_part_mut(1,2,S+1) = Matrix::from_iterator(dVdtn);
		// dn(dtn)
		for j in 2..S+1 { J.column_part_mut(j,2,S+1).copy_from_slice(dndtn[j]); }

		// dt(V) = V.(TR/P.Σ{(1-Wk/Wa).ω}+1/T.dt(T))
		// dt(T) = - Σ_ ((H-(W/Wa).Ha).w) / ( Ct.Cpa + Σ_ (Cp_Cpa)C )
		// dt(n) = Vω
		for (state, rate) in self.amounts[..specie_count-1].iter_mut().zip(production_rates.iter()) { *state = 0f64.max(*state + system.time_step * system.volume * rate); }
		let total_amount = system.pressure * system.volume / (ideal_gas_constant * self.temperature);
		self.amounts[specie_count-1] = 0f64.max(total_amount - self.amounts[..specie_count-1].iter().sum::<f64>()); // Enforces mass conservation constraint by rescaling most abundant specie (last index)
		for &amount in self.amounts.iter() { assert!(amount>= 0. && amount< total_amount,"{} {:?} {}", amount, &self.amounts, total_amount); }
		self.temperature += system.time_step * dtT;
	}
}*/
