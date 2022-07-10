use crate::math::field::Field;

pub fn get_matrix_determinant_4x4(m: &Vec<Vec<u32>>, field: &Field) -> u32 {
    let section1 = (m[0][0] * m[1][1] * m[2][2] * m[3][3]
        + m[0][0] * m[1][2] * m[2][3] * m[3][1]
        + m[0][0] * m[1][3] * m[2][1] * m[3][2])
        % field.order;

    let section2 = (m[0][0] * m[1][3] * m[2][2] * m[3][1]
        + m[0][0] * m[1][2] * m[2][1] * m[3][3]
        + m[0][0] * m[1][1] * m[2][3] * m[3][2]
        + m[0][1] * m[1][0] * m[2][2] * m[3][3]
        + m[0][2] * m[1][0] * m[2][3] * m[3][1]
        + m[0][3] * m[1][0] * m[2][1] * m[3][2])
        % field.order;

    let section3 = (m[0][3] * m[1][0] * m[2][2] * m[3][1]
        + m[0][2] * m[1][0] * m[2][1] * m[3][3]
        + m[0][1] * m[1][0] * m[2][3] * m[3][2]
        + m[0][1] * m[1][2] * m[2][0] * m[3][3]
        + m[0][2] * m[1][3] * m[2][0] * m[3][1]
        + m[0][3] * m[1][1] * m[2][0] * m[2][3])
        % field.order;

    let section4 = (m[0][3] * m[1][2] * m[2][0] * m[3][1]
        + m[0][2] * m[1][1] * m[2][0] * m[3][3]
        + m[0][1] * m[1][3] * m[2][0] * m[3][2]
        + m[0][1] * m[1][2] * m[2][3] * m[3][0]
        + m[0][2] * m[1][3] * m[2][1] * m[3][0]
        + m[0][3] * m[1][1] * m[2][2] * m[3][0])
        % field.order;

    let section5 = (m[0][3] * m[1][2] * m[2][1] * m[3][0]
        + m[0][2] * m[1][1] * m[2][3] * m[3][0]
        + m[0][1] * m[1][3] * m[2][2] * m[3][0])
        % field.order;

    field.add(
        section1,
        field.add(
            section3,
            field.subtract(section5, field.add(section2, section4)),
        ),
    )
}

pub fn get_inverse_matrix_4x4(m: &Vec<Vec<u32>>, field: &Field) -> Vec<Vec<u32>> {
    let determinant = get_matrix_determinant_4x4(m, field);
    let determinant_inverse = field.multiplicative_inverse(determinant);
    let mut adj = get_adjugate_matrix_4x4(m, field);

    for i in 0..adj.len() {
        for k in 0..adj[i].len() {
            adj[i][k] = field.multiply(adj[i][k], determinant_inverse);
        }
    }

    adj
}

pub fn matrix_multiply_4x4_1x4(m: &Vec<Vec<u32>>, v: &Vec<u32>, field: &Field) -> Vec<u32> {
    let mut res = Vec::new();

    for i in 0..v.len() {
        let mut ele = 0;
        for k in 0..m.len() {
            ele = field.add(ele, field.multiply(v[k], m[i][k]));
        }
        res.push(ele);
    }
    res
}

fn get_adjugate_matrix_4x4(m: &Vec<Vec<u32>>, field: &Field) -> Vec<Vec<u32>> {
    let mut adj: Vec<Vec<u32>> = vec![vec![0; 4]; 4];

    adj[0][0] =
        (m[1][1] * m[2][2] * m[3][3] + m[1][2] * m[2][3] * m[3][1] + m[1][3] * m[2][1] * m[3][2])
            % field.order;
    adj[0][0] = field.subtract(
        adj[0][0],
        (m[1][3] * m[2][2] * m[3][1] + m[1][2] * m[2][1] * m[3][3] + m[1][1] * m[2][3] * m[3][2])
            % field.order,
    );

    adj[0][1] =
        (m[0][3] * m[2][2] * m[3][1] + m[0][2] * m[2][1] * m[3][3] + m[0][1] * m[2][3] * m[3][2])
            % field.order;
    adj[0][1] = field.subtract(
        adj[0][1],
        (m[0][1] * m[2][2] * m[3][3] + m[0][2] * m[2][3] * m[3][1] + m[0][3] * m[2][1] * m[3][2])
            % field.order,
    );

    adj[0][2] =
        (m[0][1] * m[1][2] * m[3][3] + m[0][2] * m[1][3] * m[3][1] + m[0][3] * m[1][1] * m[3][2])
            % field.order;
    adj[0][2] = field.subtract(
        adj[0][2],
        (m[0][3] * m[1][2] * m[3][1] + m[0][2] * m[1][1] * m[3][3] + m[0][1] * m[1][3] * m[3][2])
            % field.order,
    );

    adj[0][3] =
        (m[0][3] * m[1][2] * m[2][1] + m[0][2] * m[1][1] * m[2][3] + m[0][1] * m[1][3] * m[2][2])
            % field.order;
    adj[0][3] = field.subtract(
        adj[0][3],
        (m[0][1] * m[1][2] * m[2][3] + m[0][2] * m[1][3] * m[2][1] + m[0][3] * m[1][1] * m[2][2])
            % field.order,
    );

    adj[1][0] =
        (m[1][3] * m[2][2] * m[3][0] + m[1][2] * m[2][0] * m[3][3] + m[1][0] * m[2][3] * m[3][2])
            % field.order;
    adj[1][0] = field.subtract(
        adj[1][0],
        (m[1][0] * m[2][2] * m[3][3] + m[1][2] * m[2][3] * m[3][0] + m[1][3] * m[2][0] * m[3][2])
            % field.order,
    );

    adj[1][1] =
        (m[0][0] * m[2][2] * m[3][3] + m[0][2] * m[2][3] * m[3][0] + m[0][3] * m[2][0] * m[3][2])
            % field.order;
    adj[1][1] = field.subtract(
        adj[1][1],
        (m[0][3] * m[2][2] * m[3][0] + m[0][2] * m[2][0] * m[3][3] + m[0][0] * m[2][3] * m[3][2])
            % field.order,
    );

    adj[1][2] =
        (m[0][3] * m[1][2] * m[3][0] + m[0][2] * m[1][0] * m[3][3] + m[0][0] * m[1][3] * m[3][2])
            % field.order;
    adj[1][2] = field.subtract(
        adj[1][2],
        (m[0][0] * m[1][2] * m[3][3] + m[0][2] * m[1][3] * m[3][0] + m[0][3] * m[1][0] * m[3][2])
            % field.order,
    );

    adj[1][3] =
        (m[0][0] * m[1][2] * m[2][3] + m[0][2] * m[1][3] * m[2][0] + m[0][3] * m[1][0] * m[2][2])
            % field.order;
    adj[1][3] = field.subtract(
        adj[1][3],
        (m[0][3] * m[1][2] * m[2][0] + m[0][2] * m[1][0] * m[2][3] + m[0][0] * m[1][3] * m[2][2])
            % field.order,
    );

    adj[2][0] =
        (m[1][0] * m[2][1] * m[3][3] + m[1][1] * m[2][3] * m[3][0] + m[1][3] * m[2][0] * m[3][1])
            % field.order;
    adj[2][0] = field.subtract(
        adj[2][0],
        (m[1][3] * m[2][1] * m[3][0] + m[1][1] * m[2][0] * m[3][3] + m[1][0] * m[2][3] * m[3][1])
            % field.order,
    );

    adj[2][1] =
        (m[0][3] * m[2][1] * m[3][0] + m[0][1] * m[2][0] * m[3][3] + m[0][0] * m[2][3] * m[3][1])
            % field.order;
    adj[2][1] = field.subtract(
        adj[2][1],
        (m[0][0] * m[2][1] * m[3][3] + m[0][1] * m[2][3] * m[3][0] + m[0][3] * m[2][0] * m[3][1])
            % field.order,
    );

    adj[2][2] =
        (m[0][0] * m[1][1] * m[3][3] + m[0][1] * m[1][3] * m[3][0] + m[0][3] * m[1][0] * m[3][1])
            % field.order;
    adj[2][2] = field.subtract(
        adj[2][2],
        (m[0][3] * m[1][1] * m[3][0] + m[0][1] * m[1][0] * m[3][3] + m[0][0] * m[1][3] * m[3][1])
            % field.order,
    );

    adj[2][3] =
        (m[0][3] * m[1][1] * m[2][0] + m[0][1] * m[1][0] * m[2][3] + m[0][0] * m[1][3] * m[2][1])
            % field.order;
    adj[2][3] = field.subtract(
        adj[2][3],
        (m[0][0] * m[1][1] * m[2][3] + m[0][1] * m[1][3] * m[2][0] + m[0][3] * m[1][0] * m[2][1])
            % field.order,
    );

    adj[3][0] =
        (m[1][2] * m[2][1] * m[3][0] + m[1][1] * m[2][0] * m[3][2] + m[1][0] * m[2][2] * m[3][1])
            % field.order;
    adj[3][0] = field.subtract(
        adj[3][0],
        (m[1][0] * m[2][1] * m[3][2] + m[1][1] * m[2][2] * m[3][0] + m[1][2] * m[2][0] * m[3][1])
            % field.order,
    );

    adj[3][1] =
        (m[0][0] * m[2][1] * m[3][2] + m[0][1] * m[2][2] * m[3][0] + m[0][2] * m[2][0] * m[3][1])
            % field.order;
    adj[3][1] = field.subtract(
        adj[3][1],
        (m[0][2] * m[2][1] * m[3][0] + m[0][1] * m[2][0] * m[3][2] + m[0][0] * m[2][2] * m[3][1])
            % field.order,
    );

    adj[3][2] =
        (m[0][2] * m[1][1] * m[3][0] + m[0][1] * m[1][0] * m[3][2] + m[0][0] * m[1][2] * m[3][1])
            % field.order;
    adj[3][2] = field.subtract(
        adj[3][2],
        (m[0][0] * m[1][1] * m[3][2] + m[0][1] * m[1][2] * m[3][0] + m[0][2] * m[1][0] * m[3][1])
            % field.order,
    );

    adj[3][3] =
        (m[0][0] * m[1][1] * m[2][2] + m[0][1] * m[1][2] * m[2][0] + m[0][2] * m[1][0] * m[2][1])
            % field.order;
    adj[3][3] = field.subtract(
        adj[3][3],
        (m[0][2] * m[1][1] * m[2][0] + m[0][1] * m[1][0] * m[2][2] + m[0][0] * m[1][2] * m[2][1])
            % field.order,
    );

    adj
}

#[test]
fn test_determinant() {
    let field = Field { order: 17 };

    let m: Vec<Vec<u32>> = vec![
        vec![1, 2, 6, 6],
        vec![4, 7, 3, 2],
        vec![0, 0, 0, 0],
        vec![1, 2, 2, 9],
    ];
    assert_eq!(get_matrix_determinant_4x4(&m, &field), 0);
    let m: Vec<Vec<u32>> = vec![
        vec![4, 3, 2, 2],
        vec![0, 1, 14, 3],
        vec![0, 16, 3, 3],
        vec![0, 3, 1, 1],
    ];
    assert_eq!(get_matrix_determinant_4x4(&m, &field), 15);
    let m: Vec<Vec<u32>> = vec![
        vec![1, 1, 1, 16],
        vec![1, 1, 16, 1],
        vec![1, 16, 1, 1],
        vec![16, 1, 1, 1],
    ];
    assert_eq!(get_matrix_determinant_4x4(&m, &field), 1);
}

#[test]
fn test_inverse() {
    let field = Field { order: 17 };
    let m: Vec<Vec<u32>> = vec![
        vec![1, 1, 1, 16],
        vec![1, 1, 16, 1],
        vec![1, 16, 1, 1],
        vec![16, 1, 1, 1],
    ];
    let adj = get_adjugate_matrix_4x4(&m, &field);
    let expected_adj = vec![
        vec![13, 13, 13, 4],
        vec![13, 13, 4, 13],
        vec![13, 4, 13, 13],
        vec![4, 13, 13, 13],
    ];

    for i in 0..adj.len() {
        for k in 0..adj[i].len() {
            assert_eq!(adj[i][k], expected_adj[i][k]);
        }
    }

    let m2: Vec<Vec<u32>> = vec![
        vec![1, 1, 1, 1],
        vec![1, 4, 16, 13],
        vec![1, 16, 1, 16],
        vec![1, 13, 16, 4],
    ];
    let expected_inverse: Vec<Vec<u32>> = vec![
        vec![13, 13, 13, 13],
        vec![13, 16, 4, 1],
        vec![13, 4, 13, 4],
        vec![13, 1, 4, 16],
    ];
    let actual_inverse: Vec<Vec<u32>> = get_inverse_matrix_4x4(&m2, &field);

    for i in 0..actual_inverse.len() {
        for k in 0..actual_inverse[i].len() {
            assert_eq!(actual_inverse[i][k], expected_inverse[i][k]);
        }
    }
}

#[test]
fn test_matrix_mul() {
    let field = Field { order: 17 };
    let m: Vec<Vec<u32>> = vec![
        vec![13, 13, 13, 13],
        vec![13, 16, 4, 1],
        vec![13, 4, 13, 4],
        vec![13, 1, 4, 16],
    ];
    let v: Vec<u32> = vec![3, 4, 5, 9];
    let expected_res: Vec<u32> = vec![1, 13, 3, 3];

    let actual_res = matrix_multiply_4x4_1x4(&m, &v, &field);
    for i in 0..actual_res.len() {
        assert_eq!(actual_res[i], expected_res[i]);
    }
}
