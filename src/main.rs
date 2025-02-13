use std::io;
use std::f64::consts::PI;
#[derive(Debug)]
struct Vector3D {
    x:f64,
    y:f64,
    z:f64,
}
struct PointCurvilinear{
    fi:f64,
    la:f64,
    h:f64,
}
struct PointCurvilinear1{
    h:f64,
    phi:f64,
    n:f64,
}

impl Vector3D {
    fn add(&self, other:&Vector3D)-> Vector3D{
        Vector3D {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
    fn subtract (&self, other:&Vector3D)-> Vector3D{
        Vector3D {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
    fn module(&self) -> f64 {
        ((self.x * self.x) + (self.y * self.y) + (self.z * self.z)).sqrt()
    }
}




fn read_vector(prompt: &str) -> Vector3D {
    println!("{}", prompt);
    let mut input = String::new();

    // Leer la componente x
    println!("Ingrese la componente x:");
    io::stdin().read_line(&mut input).expect("Error al leer la entrada");
    let x: f64 = input.trim().parse().expect("La componente x debe ser un número");
    input.clear(); // Limpiar el buffer de entrada

    // Leer la componente y
    println!("Ingrese la componente y:");
    io::stdin().read_line(&mut input).expect("Error al leer la entrada");
    let y: f64 = input.trim().parse().expect("La componente y debe ser un número");
    input.clear(); // Limpiar el buffer de entrada

    // Leer la componente z
    println!("Ingrese la componente z:");
    io::stdin().read_line(&mut input).expect("Error al leer la entrada");
    let z: f64 = input.trim().parse().expect("La componente z debe ser un número");

    Vector3D { x, y, z }
}
fn read_coordinate(prompt: &str) -> PointCurvilinear {
    println!("{}", prompt);
    let mut input = String::new();

    // Leer la componente x
    println!("Ingrese la componente latitud en grados decimales:");
    io::stdin().read_line(&mut input).expect("Error al leer la entrada");
    let fi: f64 = input.trim().parse().expect("La componente fi debe ser un número");
    input.clear(); // Limpiar el buffer de entrada

    // Leer la componente y
    println!("Ingrese la componente longitud en grados decimales:");
    io::stdin().read_line(&mut input).expect("Error al leer la entrada");
    let la: f64 = input.trim().parse().expect("La componente lambda debe ser un número");
    input.clear(); // Limpiar el buffer de entrada

    // Leer la componente z
    println!("Ingrese la componente altura elipsoidal en metros:");
    io::stdin().read_line(&mut input).expect("Error al leer la entrada");
    let h: f64 = input.trim().parse().expect("La componente h debe ser un número");

    PointCurvilinear { fi, la, h }
}
fn compute_n(a: f64, b: f64, p_geo: &PointCurvilinear1) -> f64 {
    let e2 = 1.0 - (b.powi(2) / a.powi(2));
    a / (1.0 - e2 * p_geo.phi.sin().powi(2)).sqrt()
}
fn computeN (a:f64, b:f64, punto:&PointCurvilinear)->f64{
    let fi_radianes = punto.fi * (PI / 180.0); // Convertir grados a radianes
    let numerador = a * a;
    println!("numerador es: {}", numerador);
    
    let denominador = (((a * a) * fi_radianes.cos().powi(2)) + ((b * b) * fi_radianes.sin().powi(2))).sqrt();
    println!("denominador es: {}", denominador);
    numerador / denominador
}
fn computeCartesian(punto:&PointCurvilinear,primerVertical:f64, a:f64, b:f64)->Vector3D{
    let fi_radianes = punto.fi * (PI / 180.0);
    let lambda_radianes = punto.la * (PI / 180.0);
    let x = (primerVertical+punto.h)* fi_radianes.cos()*lambda_radianes.cos();
    let y = (primerVertical+punto.h)* fi_radianes.cos()*lambda_radianes.sin();
    let z = ((primerVertical * b.powi(2) / a.powi(2))+punto.h)* fi_radianes.sin();
    Vector3D {x,y,z}
}
fn computeLongitude(punto: &Vector3D)-> f64{
    let y = punto.y;
    let x = punto.x;
    let arctan2_yx = y.atan2(x);
    arctan2_yx 
}

fn compute_iterative(punto: &Vector3D, a: f64, b: f64) -> (f64, f64, f64) {
    // Valores iniciales
    let mut n = a;
    let z = punto.z;
    let x = punto.x;
    let y = punto.y;
    let e2 = 1.0 - (b.powi(2) / a.powi(2));
    let p = (x.powi(2) + y.powi(2)).sqrt();
    let magnitud = punto.module();
    let mut h = magnitud - (a * b).sqrt();
    let phi1 = z / p;
    let mut phi2 = 1.0 - (e2 * n / (n + h));
    let mut phi3 = phi1 / phi2;
    let mut phi = phi3.atan();
    let mut contador =0;

    // Umbral de convergencia
    let threshold = 1e-6;

    // Bucle iterativo
    loop {
        contador = contador + 1;
        // Guardar los valores anteriores
        let h_prev = h;
        let phi_prev = phi;

        // Calcular nuevos valores
        n = compute_n(a, b, &PointCurvilinear1 { h, phi, n });
        h = compute_h_iterative(punto, &PointCurvilinear1 { h, phi, n }, phi, a, b);
        phi = compute_latitude_iterative(punto, &PointCurvilinear1 { h, phi, n }, h, a, b);
        println!("iteracion: {} ", contador);
        println!("phi: {} grados", phi * 180.0 / PI);
        println!("h: {} ", h);
        println!("N: {} ", n);
        // Verificar convergencia
        if (h - h_prev).abs() < threshold && (phi - phi_prev).abs() < threshold {
            break;
        }
    }

    (h, phi, n)
}

// Función para calcular h iterativamente
fn compute_h_iterative(punto: &Vector3D, p_geo: &PointCurvilinear1, phi: f64, a: f64, b: f64) -> f64 {
    let z = punto.z;
    let x = punto.x;
    let y = punto.y;
    let fphi = phi; // phi ya está en radianes
    let p = (x.powi(2) + y.powi(2)).sqrt();
    let n = compute_n(a, b, p_geo);
    let h1 = p / fphi.cos();
    h1 - n
}

// Función para calcular la latitud iterativamente
fn compute_latitude_iterative(punto: &Vector3D, p_geo: &PointCurvilinear1, h: f64, a: f64, b: f64) -> f64 {
    let z = punto.z;
    let x = punto.x;
    let y = punto.y;
    let e2 = 1.0 - (b.powi(2) / a.powi(2));
    let p = (x.powi(2) + y.powi(2)).sqrt();
    let n = compute_n(a, b, p_geo);
    let z_p = z / p;
    let temp_e = 1.0 - (e2 * n / (n + h));
    let temp_inv = 1.0 / temp_e;
    let ze = z_p * temp_inv;
    ze.atan()
}


fn main() {
    //WGS84 parameters
    let a = 6378137.000;
    let b = 6356752.3142;
    // Pedir al usuario los dos vectores
    //let v1 = read_vector("Ingrese el primer vector:");
    //let v2 = read_vector("Ingrese el segundo vector:");
    //Pedir al usuario dos coordenadas 
    let p1  =read_coordinate("Ingrese las coordenadas del primer punto");
    let p2  =read_coordinate("Ingrese las coordenadas del segundo punto");
    //Calcular los vectores Normal de cada punto 
    let n1 = computeN(a,b,&p1);
    let n2 = computeN(a,b,&p2);
    //Obtener las coordenadas cartesianas de cada punto
    let v1 = computeCartesian(&p1,n1,a,b);
    let v2 = computeCartesian(&p2,n2,a,b);
    let resultado_resta = v1.subtract(&v2);
    let modulo = resultado_resta.module();
    println!("La primer vertical n1 es: {}", n1);
    println!("La segunda vertical n1 es: {}", n2);
    println!("Las coordenadas cartesianas del primer punto son: {:?}", v1);
    println!("Las coordenadas cartesianas del segundo punto son: {:?}", v2);
    println!("El vector formado de ambos puntos es: {:?}", resultado_resta);
    println!("El módulo del vector resultante es: {}", modulo);
    println!("Calculando las coordenadas geodesicas de vuelta");
    let longitud1 = computeLongitude(&v1);
    println!("longitud1: {}", longitud1);
    
    let longitud2 =computeLongitude(&v2);
    println!("longitud2: {}", longitud2);
    let (h11, phi11, n11) = compute_iterative(&v1, a, b);
    println!("h: {}", h11);
    println!("phi: {} radianes", phi11);
    println!("phi: {} grados", phi11 * 180.0 / PI);
    println!("n: {}", n11);
    let (h22, phi22, n22) = compute_iterative(&v2, a, b);
    println!("h: {}", h22);
    println!("phi: {} radianes", phi22);
    println!("phi: {} grados", phi22 * 180.0 / PI);
    println!("n: {}", n22);



}