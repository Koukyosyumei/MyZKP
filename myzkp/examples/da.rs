use myzkp::modules::das::{
    avail::Avail,
    celestia::Celestia,
    eigenda::EigenDA,
    utils::{reset_metrics, DataAvailabilitySystem, SamplePosition, METRICS},
};

fn main() {
    let data_size = 64;
    let sqrt_data_size = 8;
    let num_sampling = 26;

    let data: Vec<_> = (0..63).collect();

    println!("# EigenDA");
    let params = EigenDA::setup(16, 4.0);
    let encoded = EigenDA::encode(&data, &params);
    let commit = EigenDA::commit(&encoded, &params);
    for i in 0..9 {
        let position = SamplePosition {
            row: 0,
            col: i,
            is_row: false,
        };
        let isvalid = EigenDA::verify(&position, &encoded, &commit, &params);
        assert!(isvalid);
    }
    METRICS.with(|metrics| {
        println!("{:#?}", *metrics.borrow());
    });
    reset_metrics();

    println!("# Celestia");
    let params = Celestia::setup(8, 2.0);
    let encoded = Celestia::encode(&data, &params);
    let commit = Celestia::commit(&encoded, &params);
    for i in 0..25 {
        let position = SamplePosition {
            row: i / 16,
            col: i % 16,
            is_row: false,
        };
        assert!(Celestia::verify(&position, &encoded, &commit, &params));
    }
    METRICS.with(|metrics| {
        println!("{:#?}", *metrics.borrow());
    });
    reset_metrics();

    println!("# Avail");
    let params = Avail::setup(8, 2.0);
    let encoded = Avail::encode(&data, &params);
    let commit = Avail::commit(&encoded, &params);
    for i in 0..25 {
        let position0 = SamplePosition {
            row: i / 16,
            col: i % 16,
            is_row: false,
        };
        assert!(Avail::verify(&position0, &encoded, &commit, &params));
    }
    METRICS.with(|metrics| {
        println!("{:#?}", *metrics.borrow());
    });
}
