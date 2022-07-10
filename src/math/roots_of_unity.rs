use crate::math::field::Field;
pub fn get_roots_of_unity(n: u32, field: &Field) -> Vec<u32> {
    let mut res = Vec::new();
    for i in 0..field.order {
        if field.exponent(i, n) == 1 {
            res.push(i);
        }
    }

    // Manually switch last two to be in line with plonk_by_hand article
    res.swap(3, 2);

    res
}

pub fn get_coset(k: u32, roots: &Vec<u32>, field: &Field) -> Vec<u32> {
    roots.iter().map(|root| field.multiply(*root, k)).collect()
}

#[test]
fn test_roots() {
    let expected_roots = vec![1, 4, 16, 13];
    let actual_roots = get_roots_of_unity(4, &Field { order: 17 });
    assert_eq!(expected_roots, actual_roots);
}

#[test]
fn test_cosets() {
    let field = Field { order: 17 };
    let expected_coset_1 = vec![2, 8, 15, 9];
    let actual_coset_1 = get_coset(2, &vec![1, 4, 16, 13], &field);
    assert_eq!(expected_coset_1, actual_coset_1);

    let expected_coset_2 = vec![3, 12, 14, 5];
    let actual_coset_2 = get_coset(3, &vec![1, 4, 16, 13], &field);
    assert_eq!(expected_coset_2, actual_coset_2);
}
