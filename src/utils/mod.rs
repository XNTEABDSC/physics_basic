
#[macro_export]
macro_rules! derive_zero {
	($t:tt) => {
impl <Num: num_traits::identities::Zero> num_traits::identities::Zero for $t<Num> {
	fn zero() -> Self {
		Self(Num::zero())
	}

	fn is_zero(&self) -> bool {
		Num::is_zero(&self.0)
	}
}
impl <Num: num_traits::identities::Zero> core::default::Default for $t<Num> {
    fn default() -> Self {
		Self(Num::zero())
	}
}
	};
}