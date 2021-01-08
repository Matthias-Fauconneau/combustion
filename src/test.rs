use super::*;

#[test] fn test() {
    let system = std::fs::read("CH4+O2.ron").unwrap();
    let Simulation{species_names, system, ..} = Simulation::<35>::new(&system).unwrap();
    let transport = |single_specie, temperature_C| {
        let pressure_R = 1e5/(kB*NA);
        let temperature = 273.15+temperature_C;
        let amount = pressure_R / temperature * System::<35>::volume;
        system.transport(pressure_R, &State{temperature, amounts: eval(species_names, |specie| if specie==single_specie {amount} else {0.})})
    };
    let viscosity = |single_specie, T, expected| { let e = f64::abs(transport(single_specie, T).viscosity*1e6-expected)/expected; println!("{}", e); assert!(e < 0.07); };
    viscosity("Ar",25., 22.58);
    viscosity("CH4",25., 11.07);
    let thermal_conductivity = |single_specie, T, expected| { let got = transport(single_specie, T).thermal_conductivity; let e = f64::abs(got-expected)/expected; println!("{}", e);  assert!(e < 0.05, "{} {}", got, expected); };
    thermal_conductivity("Ar",500., 0.03650);
}
