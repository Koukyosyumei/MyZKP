use myzkp::modules::das::{
    celestia::Celestia,
    utils::{DataAvailabilitySystem, SamplePosition, METRICS},
};

fn main() {
    let params = Celestia::setup(8, 2.0);

    let data: Vec<_> = (0..63).collect();
    let encoded = Celestia::encode(&data, &params);
    let commit = Celestia::commit(&encoded, &params);

    for i in 0..10 {
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
}
