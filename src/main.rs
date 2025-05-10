use std::io;
use std::f64::consts::PI;
use serde::Deserialize;
use std::fs::File;
use std::error::Error;
use std::path::Path;
//#[derive(Debug)]
#[derive(Debug, Clone, Copy)]
struct Vector3D {
    x:f64,
    y:f64,
    z:f64,
}
struct Datadist {
    azimuth12:f64,
    dist12:f64,
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
//Estructura axuliar para deserializar lo que leemos del archivo que estamos leyendo
#[derive(Debug, Deserialize)]
struct CsvPoint{
    #[serde(rename = "Longitude")]
    longitud:f64,
    #[serde(rename = "Latitude")]
    latitud:f64,
    #[serde(rename ="Ellipsoidal height")]
    elipsoidal_height: f64,
}

pub fn leer_puntos_archivo<P: AsRef<Path>>(file_path: P) -> Result<Vec<PointCurvilinear>, Box<dyn Error>> {
    let file = File::open(file_path)?;
    let mut rdr = csv::Reader::from_reader(file);
    let mut points = Vec::new();

    for result in rdr.deserialize() {
        let record: CsvPoint = result?;
        let point = PointCurvilinear {
            fi: record.latitud,
            la: record.longitud,
            h: record.elipsoidal_height,
        };
        points.push(point);
    }

    // Verificación opcional: asegurar que se leyó al menos un punto
    if points.is_empty() {
        eprintln!("Advertencia: El archivo CSV no contenía puntos válidos");
    }

    Ok(points)
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


fn read_azdist(prompt: &str) -> Datadist {
    println!("{}", prompt);
    let mut input = String::new();
    println!("Ingrese la distancia al punto 2:");
    io::stdin().read_line(&mut input).expect("Error al leer la entrada");
    let dist12: f64 = input.trim().parse().expect("La distancia debe ser un número");
    input.clear(); // Limpiar el buffer de entrada
    println!("Ingrese el azimuth al punto 2:");
    io::stdin().read_line(&mut input).expect("Error al leer la entrada");
    let azimuth12: f64 = input.trim().parse().expect("el azimuth debe ser un numero");
    input.clear(); // Limpiar el buffer de entrada
    Datadist{ azimuth12,dist12 }
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
fn compute_nsencillo(a: f64, b: f64, latitud: f64) -> f64 {
    let e2 = 1.0 - (b.powi(2) / a.powi(2));
    a / (1.0 - e2 * latitud.sin().powi(2)).sqrt()
}
fn compute_msencillo(a: f64, b: f64, latitud:f64) -> f64 {
    let e2 = 1.0 - (b.powi(2) / a.powi(2));
    (a*(1.0 - e2)) / (1.0 - e2 * latitud.sin().powi(2)).powf(1.5)
}
fn compute_m(a: f64, b: f64, p_geo: &PointCurvilinear) -> f64 {
    let e2 = 1.0 - (b.powi(2) / a.powi(2));
    (a*(1.0 - e2)) / (1.0 - e2 * p_geo.fi.sin().powi(2)).powf(1.5)
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
//Funcion para calcular el arco de meridiano entre 2 latitudes diferentes usando aproximacion por expansion de Binomio de Newton y solucionando de manera analitica
fn arcoMeridiano(p_geo1: &PointCurvilinear, p_geo2: &PointCurvilinear, a: f64, b:f64)->f64{

    let phi1_rad = p_geo1.fi* (PI / 180.0);
    let phi2_rad = p_geo2.fi* (PI / 180.0);
    let delta_phi = phi2_rad - phi1_rad;
    let e2 = 1.0 - (b.powi(2) / a.powi(2));

    // Términos de la expansión en serie
    let term1 = delta_phi;
    let term2 = (3.0 / 4.0) * e2 * delta_phi;
    let term3 = -(3.0 / 8.0) * e2 * (f64::sin(2.0 * phi2_rad) - f64::sin(2.0 * phi1_rad));
    let term4 = (15.0 / 64.0) * e2.powi(2) * delta_phi;
    let term5 = -(15.0 / 32.0) * e2.powi(2) * (f64::sin(2.0 * phi2_rad) - f64::sin(2.0 * phi1_rad));
    let term6 = (15.0 / 256.0) * e2.powi(2) * (f64::sin(4.0 * phi2_rad) - f64::sin(4.0 * phi1_rad));

    // Sumar todos los términos
    let s = a * (1.0 - e2) * (term1 + term2 + term3 + term4 + term5 + term6);
    s
}
fn arcoMeridianoPromedio(p_geo1: &PointCurvilinear, p_geo2: &PointCurvilinear, a: f64, b:f64)->f64{

    let phi1_rad = p_geo1.fi* (PI / 180.0);
    let phi2_rad = p_geo2.fi* (PI / 180.0);
    let phimedia = (phi1_rad+phi2_rad)/2.0;
    let delta_phi = phi2_rad - phi1_rad;
    let e2 = 1.0 - (b.powi(2) / a.powi(2));
    let m = compute_msencillo(a,b,phimedia);
    let s= delta_phi*m;
    s


   
}

//Funcion para calcular el arco de paralelo ente 2 longitudes diferentes radio de circulo utilizando la primer vertical o N
fn arcoParalelo(p_geo1: &PointCurvilinear, p_geo2: &PointCurvilinear,a: f64, b:f64)->f64{
    let phi1_rad = p_geo1.fi* (PI / 180.0);
    let phi2_rad = p_geo2.fi* (PI / 180.0);
    let la1_rad = p_geo1.la * (PI / 180.0);
    let la2_rad = p_geo2.la * (PI / 180.0);
    let menor = phi1_rad.min(phi2_rad);
    //Utilizo la latiud menor para sacar el radio de la primer normal porque es el cateto adyacente del triangulo rectangulo formado entre ambos puntos
    let nm = compute_nsencillo(a,b,menor);
    println!("nm {}", nm);
    let radio = nm*menor.cos();
    //diferencia de longitudes
    let deltaLambda=  la2_rad-la1_rad;
    let s = radio * deltaLambda;
    s
}

//Funcio para calcular distancia entre dos puntos usando el arco de paralelo y el arco de meridiano como los catetos de un  triangulo rectangulo
fn calcularDistanciaArcos(arcoM: f64, arcoP: f64)-> f64{
    let s = arcoM.powi(2)+ arcoP.powi(2);
    let s2 =s.sqrt();
    s2
}
// Funcion para calcular distancia y azimuth con Puissant a traves de iteraciones
fn computePuissantInverso(p_geo1: &PointCurvilinear,p_geo2: &PointCurvilinear, a: f64, b: f64)-> (f64, f64) {
    let n1 = computeN(a, b, p_geo1);
    let m1 = compute_m(a, b, p_geo1);
    let n2 = computeN(a, b, p_geo2);
    let phi1=p_geo1.fi* (PI / 180.0);
    let phi2=p_geo2.fi* (PI / 180.0);
    let lam1=p_geo1.la* (PI / 180.0);
    let lam2=p_geo2.la* (PI / 180.0);
    let delta_lambda = (lam1 - lam2).abs();
    let delta_phi = (phi1 - phi2).abs();
    let e2 = 1.0 - (b.powi(2) / a.powi(2));
    // primera iteracion, se usa primer termino para aproximar el azimuth directo
    let primeralpha1216= delta_phi * m1 / (1.0 - (3.0 * e2 * phi1.sin() * phi1.cos()) / (2.0 * (1.0 - e2 * phi1.sin().powi(2))));    
    let complemento_phi = (PI/2.0)-phi2;
    let primeralpha1217 =delta_lambda*n2*complemento_phi.sin();
    let mut cocientePrimero = (primeralpha1217 / primeralpha1216).atan();
    // primera distancia
    let mut primers1217 = delta_lambda*n2*complemento_phi.sin()/cocientePrimero;
    println!("Azimuth primera aprox : {}", cocientePrimero* (180.0 / PI));
    println!("Distancia primera aprox : {}", primers1217);  
    // iteraciones
    let thresholdAzimuth = 1e-6;
    let thresholdDistancia = 1e-3;
    let mut iter_count = 0;
    let max_iter = 15;
    loop {
        iter_count += 1;
        if iter_count > max_iter {
            println!("Advertencia: No converge después de {} iteraciones", max_iter);
            break;
        }
        let term2 = primers1217.powi(2) * phi1.tan() * primeralpha1217.sin().powi(2) / (2.0 * m1 * n1);
        let term3 = primers1217.powi(3) * primeralpha1217.cos() * primeralpha1217.sin().powi(2) * (1.0 + 3.0 * phi1.tan().powi(2)) / (6.0 * m1 * n1.powi(2));
        let alpha1216 = ((delta_phi / (1.0 - (3.0 * e2 * phi1.sin() * phi1.cos()) / (2.0 * (1.0 - e2 * phi1.sin().powi(2)))))+ term2 + term3) * m1;
        let term2l = (primers1217.powi(3)/(6.0*n2.powi(3))) * (primeralpha1217.sin()/complemento_phi.sin());
        let term3l = (primers1217.powi(3)/(6.0*n2.powi(3))) * (primeralpha1217.sin().powi(3)/complemento_phi.sin().powi(3));
        let alpha1217 = (delta_lambda + term2l - term3l)*n2*complemento_phi.sin();
        let cocientesegundo = (alpha1217 / alpha1216).atan();
        let segundos1217  = (delta_lambda + term2l - term3l)*n2 /  (cocientesegundo.sin()/complemento_phi.sin());
        println!("Iteración: Azimuth: {}, Distancia: {}", cocientesegundo * (180.0 / PI), segundos1217);

        if (cocientesegundo - cocientePrimero).abs() < thresholdAzimuth && (segundos1217 - primers1217).abs()  < thresholdDistancia{
            cocientePrimero = cocientesegundo;
            primers1217 =segundos1217;
            break;
        }
        cocientePrimero = cocientesegundo;
        primers1217 =segundos1217;   

    }

    (cocientePrimero,primers1217)


}

// Funcion para calcular el punto siguiente con la formula de Puissant hasta 150 kilometros
fn compute_p2_Puissant(p_geo: &PointCurvilinear, a: f64, b: f64, s: f64, azimuth: f64)-> PointCurvilinear {
    let n = computeN(a, b, p_geo);
    let m = compute_m(a, b, p_geo);
    let e2 = 1.0 - (b.powi(2) / a.powi(2));
    let phi1=p_geo.fi* (PI / 180.0);
    let lam1=p_geo.la* (PI / 180.0);
    let s12=s;
    let alpha12 = azimuth* (PI / 180.0);
    //Primero calculamos delta phi para obtener phi2
    let term1 = s12 * alpha12.cos() / m;
    let term2 = s12.powi(2) * phi1.tan() * alpha12.sin().powi(2) / (2.0 * m * n);
    let term3 = s12.powi(3) * alpha12.cos() * alpha12.sin().powi(2) * (1.0 + 3.0 * phi1.tan().powi(2)) / (6.0 * m * n.powi(2));
    let rho_prime_prime = term1 - term2 - term3;
    //let s1216 = ((delta_phi / (1.0 - (3.0 * e2 * phi1.sin() * phi1.cos()) / (2.0 * (1.0 - e2 * phi1.sin().powi(2))))) + term2 + term3)* m/ alpha12.cos();
    
    let delta_phi = rho_prime_prime * (1.0 - (3.0 * e2 * phi1.sin() * phi1.cos()) / (2.0 * (1.0 - e2 * phi1.sin().powi(2))));
    //let alpha1216 = ((delta_phi / (1.0 - (3.0 * e2 * phi1.sin() * phi1.cos()) / (2.0 * (1.0 - e2 * phi1.sin().powi(2)))))+ term2 + term3) * m;
    //let primeralpha1216= delta_phi * m / (1.0 - (3.0 * e2 * phi1.sin() * phi1.cos()) / (2.0 * (1.0 - e2 * phi1.sin().powi(2))));
    let phi2 = delta_phi* (180.0 / PI) + p_geo.fi;
    println!("La delta Phi es de : {}", delta_phi* (180.0 / PI));
    //calculamos delta lambda para obtener lambda2
    let latitud_media = (p_geo.fi+phi2)/2.0;
    let nm = compute_nsencillo(a,b,phi2);
    let complemento_phi = (PI/2.0)-phi1;
    let term1l = (s12/nm) * (alpha12.sin()/complemento_phi.sin());
    let term2l = (s12.powi(3)/(6.0*nm.powi(3))) * (alpha12.sin()/complemento_phi.sin());
    let term3l = (s12.powi(3)/(6.0*nm.powi(3))) * (alpha12.sin().powi(3)/complemento_phi.sin().powi(3));
    
    //let s1217  = (delta_lambda + term2l - term3l)*nm /  (alpha12.sin()/complemento_phi.sin());
    //let alpha1217 = (delta_lambda + term2l - term3l)*nm /  (1/complemento_phi.sin());
    //let primeralpha1217 =delta_lambda*nm*complemento_phi.sin();
    //let primers1217 = delta_lambda*nm*complemento_phi.sin()/primeralpha1217;


    let delta_lambda = term1l-term2l+term3l;
    println!("La delta Lambda es de : {}", delta_lambda* (180.0 / PI));
    let lambda2 = delta_lambda* (180.0 / PI) + p_geo.la;
    let height =0.0;
    PointCurvilinear{fi:phi2,la:lambda2,h:p_geo.h}
}

fn latitudB1(p_geo: &PointCurvilinear,a: f64, b: f64)-> f64{
    let phi1=p_geo.fi* (PI / 180.0);
    let e2 = 1.0 - (b.powi(2) / a.powi(2));
    let b1 =((1.0_f64-e2).sqrt() * phi1.tan()).atan();
    b1
}

/*
fn lonBessel(
    a: f64,
    b: f64,
    sigmat: f64,
    azimuth12: f64,
    b0: f64,
    b1: f64,
    b2: f64
) -> f64 {
    let threshold = 0.00001;
    let e2 = 1.0 - (b.powi(2) / a.powi(2));

    // Paso 1: calcular lambda inicial (esférico)
    let mut lambda_prev = ((sigmat.sin() * azimuth12.sin()) / b2.cos()).asin();

    // Variables de control
   
    // Azimut medio (alpha_m)
    let alpha = ((b1.cos() * b2.cos() * lambda_prev.sin()) / sigmat.sin()).asin();

    let c2sigmam = sigmat.cos() - ((2.0 * b1.sin() * b2.sin()) / alpha.cos().powi(2));
    let c4sigmam = 2.0 * c2sigmam.powi(2) - 1.0;

    // Coeficientes de serie
    let aprima = 1.0
        + (e2 / 4.0)
        + (e2.powi(2) / 8.0)
        - (e2 * b0.sin().powi(2) / 8.0)
        - (e2.powi(2) * b0.sin().powi(2) / 8.0)
        + (3.0 * e2.powi(2) * b0.sin().powi(4) / 64.0);

    let bprima = (e2 * b0.sin().powi(2) / 8.0)
        + (e2.powi(2) * b0.sin().powi(2) / 8.0)
        - (e2.powi(2) * b0.sin().powi(4) / 16.0);

    let cprima = e2.powi(2) * b0.sin().powi(4) / 64.0;

    // Calcular L nuevo
    let mut l_prev = (e2 / 2.0)
        * b0.cos()
        * (aprima * sigmat
            + bprima * sigmat.sin() * c2sigmam
            + cprima * (2.0 * sigmat).sin() * c4sigmam / 2.0);


    let mut eleMayuscula_prev = lambda_prev - l_prev;
    let mut dif_prev = lambda_prev-eleMayuscula_prev; //primera iteracion
    println!("dif_prev es: {}", dif_prev);
    let mut iter_count = 0;
    let max_iter = 15;

    loop {
        iter_count += 1;
        if iter_count > max_iter {
            println!("Advertencia: No converge después de {} iteraciones", max_iter);
            break;
        }
        let sigmai = (b1.sin() * b2.sin() + (b1.cos() * b2.cos() * dif_prev)).acos();
        println!("sigmai es: {}", sigmai);
        let lambda_new = ((sigmai.sin() * azimuth12.sin()) / b2.cos()).asin();
        println!("lambda_new es: {}", lambda_new);
        // Azimut medio (alpha_m)
        println!("sigmat.sin() es: {}", sigmat.sin());
        println!("b1.cos()es: {}", b1.cos());
        println!("b3.cos()es: {}", b2.cos());
        println!("lambda new sin() es: {}", lambda_new.sin());
        let alpha = ((b1.cos() * b2.cos() * lambda_new.sin()) / sigmai.sin()).asin();
        println!("alpha es: {}", alpha);

        let c2sigmam = sigmai.cos() - ((2.0 * b1.sin() * b2.sin()) / alpha.cos().powi(2));
        let c4sigmam = 2.0 * c2sigmam.powi(2) - 1.0;

        // Coeficientes de serie
        let aprima = 1.0
            + (e2 / 4.0)
            + (e2.powi(2) / 8.0)
            - (e2 * b0.sin().powi(2) / 8.0)
            - (e2.powi(2) * b0.sin().powi(2) / 8.0)
            + (3.0 * e2.powi(2) * b0.sin().powi(4) / 64.0);

        let bprima = (e2 * b0.sin().powi(2) / 8.0)
            + (e2.powi(2) * b0.sin().powi(2) / 8.0)
            - (e2.powi(2) * b0.sin().powi(4) / 16.0);

        let cprima = e2.powi(2) * b0.sin().powi(4) / 64.0;

        // Calcular L nuevo
        let l_new = (e2 / 2.0)
            * b0.cos()
            * (aprima * sigmai
                + bprima * sigmai.sin() * c2sigmam
                + cprima * (2.0 * sigmai).sin() * c4sigmam / 2.0);
        println!("l_new es: {}", l_new);


        let eleMayuscula_new = lambda_new - l_new;
        let dif_new = lambda_new-eleMayuscula_new;
        // Verificar convergencia de (lambda - L)
        println!("dif_new es: {}", dif_new);
        println!("dif_prev es: {}", dif_prev);
        let difdifs = (dif_prev - dif_new).abs();
        println!("difs es: {}", difdifs);
        if (dif_prev - dif_new).abs()  < threshold {
            eleMayuscula_prev= eleMayuscula_new;
            l_prev = l_new;
            break;
        }
        dif_prev = dif_new;
        l_prev = l_new;
        eleMayuscula_prev= eleMayuscula_new;
    }
    println!("eleMayuscula_prev es: {}", eleMayuscula_prev);
    eleMayuscula_prev  // Esta es la diferencia de longitud elipsoidal L
}*/

fn lonBessel(
    a: f64,
    b: f64,
    sigmat: f64,
    azimuth12: f64,
    b0: f64,
    b1: f64,
    b2: f64,
) -> (f64,f64) {
    let threshold = 1e-12; // Tolerancia más estricta para precisión geodésica
    let e2 = (a.powi(2) - b.powi(2)) / a.powi(2); // Excentricidad al cuadrado correcta

    // Coeficientes constantes (calculados una sola vez)
    let aprima = 1.0 
        + e2/4.0 
        + e2.powi(2)/8.0 
        - (e2/8.0 + e2.powi(2)/8.0)*b0.sin().powi(2) 
        + (3.0*e2.powi(2)/64.0)*b0.sin().powi(4);

    let bprima = (e2/8.0 + e2.powi(2)/8.0)*b0.sin().powi(2) 
        - (e2.powi(2)/16.0)*b0.sin().powi(4);

    let cprima = (e2.powi(2)/64.0)*b0.sin().powi(4);

    // Paso 1: Calcular lambda inicial (aproximación esférica)
    let mut lambda = ((sigmat.sin() * azimuth12.sin()) / b2.cos()).asin();
    let mut l_diff;
    let mut l_result = 0.0;
    let mut alpha_result =0.0;

    // Iteración principal
    for _ in 0..15 { // Límite de iteraciones razonable
        // Calcular sigma_i usando el lambda actual
        let cos_sigmai = b1.sin()*b2.sin() + b1.cos()*b2.cos()*lambda.cos();
        let sigmai = cos_sigmai.acos();

        // Calcular términos trigonométricos medios
        let alpha = ((b1.cos() * b2.cos() * lambda.sin()) / sigmai.sin()).asin();
        alpha_result= alpha;
        let cos_2sigma_m = sigmai.cos() - (2.0*b1.sin()*b2.sin())/alpha.cos().powi(2);
        let cos_4sigma_m = 2.0*cos_2sigma_m.powi(2) - 1.0;

        // Calcular (λ - L) según ecuación (179)
        l_diff = (e2/2.0) * b0.cos() * (
            aprima * sigmai
            + bprima * sigmai.sin() * cos_2sigma_m
            + (cprima/2.0) * (2.0*sigmai).sin() * cos_4sigma_m
        );

        // Actualizar L = λ - (λ - L)
        l_result = lambda - l_diff;

        // Verificar convergencia
        if (lambda - l_result).abs() < threshold {
            break;
        }

        // Preparar siguiente iteración
        lambda = l_result;
    }

    (l_result,alpha_result)
}

fn iterarSigma(a:f64,b:f64, s:f64, b1:f64, azimuth:f64)->f64{
    //Primera iteracion para encotnrar sigma, primer termino de ecuacion 166
    let sigma1 = (b1.tan() / azimuth.cos()).atan(); //sigma1 con eeucacion 144a
    println!("Sigma 1  es: {}", sigma1);
    let b0 = (azimuth.sin()*b1.cos()).acos();//b0 latitud reducida en el punto de retorno
    //Calculo de las constantes de la ecucion 165
    let ePrima2 =(a.powi(2)-b.powi(2)) / b.powi(2);
    println!("eprima2  es: {}", ePrima2);
    let k2 = ePrima2 * b0.sin().powi(2);
    println!("K  es: {}", k2);
    let aMayor = 1.0+(k2/4.0)-(3.0*(k2.powi(2))/64.0);
    println!("A  es: {}", aMayor);
    let bMayor = (k2/4.0)-(k2.powi(2)/16.0);
    println!("B  es: {}", bMayor);
    let cMayor = (k2.powi(2)/64.0);
    println!("C  es: {}", cMayor);    
    //principio de las iteraciones sigma0 es igual al primer termino de la ecuacion solamente
    let mut sigmai= s/(aMayor*b);
    let thresholdSigma=0.00001;
    //let mut dosSigmam = 2*sigma1+sigmai;
    //iteraciones
    loop {
        println!("Sigma i  es: {}", sigmai);
        let dosSigmam = 2.0*sigma1+sigmai;
        println!("2sisgma m  es: {}", dosSigmam);
        let term1 = s/(aMayor*b);
        let term2 = (bMayor*((dosSigmam.cos())*(sigmai.sin())))/aMayor;
        println!("termino 2  es: {}", term2);
        let term3 = (cMayor*(((2.0*dosSigmam).cos())*((2.0*sigmai).sin()))/(2.0*aMayor));
        println!("termino 3  es: {}", term3);
        let sigmaimas1 = term1+ term2+term3;
        if (sigmaimas1 - sigmai).abs()<thresholdSigma{
            sigmai= sigmaimas1;
            break;
        }
        sigmai= sigmaimas1;
    }
    sigmai
    }

fn reverse_azimuthBessel(azimuth:f64, b2:f64)->f64{
    let cociente = azimuth.sin() / b2.cos();
    let azimuthreverse = cociente.asin();
    azimuthreverse
}

fn compute_bessel_inverse(
    p_geo1: &PointCurvilinear,
    p_geo2: &PointCurvilinear,
    a: f64,
    b: f64
) -> (f64, f64) {
    // 1. Conversión a radianes y parámetros iniciales
    let phi1 = p_geo1.fi.to_radians();
    let phi2 = p_geo2.fi.to_radians();
    let lam1 = p_geo1.la.to_radians();
    let lam2 = p_geo2.la.to_radians();
    
    let e2 = (a.powi(2) - b.powi(2)) / a.powi(2); // Excentricidad al cuadrado
    let e_prime2 = (a.powi(2) - b.powi(2)) / b.powi(2); // Segunda excentricidad
    
    // 2. Cálculo de latitudes reducidas (β)
    let beta1 = ((1.0 - e2).sqrt() * phi1.tan()).atan();
    let beta2 = ((1.0 - e2).sqrt() * phi2.tan()).atan();
    
    // 3. Diferencia de longitud y cálculo inicial de σ
    let delta_lambda = (lam2 - lam1).abs();
    let cos_sigma = beta1.sin() * beta2.sin() + beta1.cos() * beta2.cos() * delta_lambda.cos();
    let mut sigma = cos_sigma.acos();
    
    // 4. Cálculo inicial del azimut α
    let alpha = ((beta1.cos() * beta2.cos() * delta_lambda.sin()) / sigma.sin()).asin();
    
    // 5. Coeficientes constantes (calculados una sola vez)
    let b0 = (alpha.sin() * beta1.cos()).acos();
    let k2 = e_prime2 * alpha.cos().powi(2);
    
    let a_coeff = 1.0 + (k2/4.0) - (3.0*k2.powi(2)/64.0);
    let b_coeff = (k2/4.0) - (k2.powi(2)/16.0);
    let c_coeff = (k2.powi(2)/64.0);
    
    // 6. Iteración para convergencia (ecuación 179 del documento)
    let mut iter_count = 0;
    let max_iter = 15;
    let threshold = 1e-12;
    let mut lambda_diff_prev;
    let mut lambda_prev = delta_lambda;
    let mut sigma_i = 0.0;
    loop {
        iter_count += 1;
        
        // 7. Cálculo de términos intermedios
        let cos_sigma_i = beta1.sin() * beta2.sin() + beta1.cos() * beta2.cos() * lambda_prev.cos();
        sigma_i = cos_sigma_i.acos();
        
        let alpha_i = ((beta1.cos() * beta2.cos() * lambda_prev.sin()) / sigma_i.sin()).asin();
        let b0_i = (alpha_i.sin() * beta1.cos()).acos();
        
        let two_sigma_m = 2.0 * ((beta1.sin() * beta2.sin() / sigma_i.cos()).atan()) - sigma_i;
        let cos_2sigma_m = two_sigma_m.cos();
        let cos_4sigma_m = 2.0 * cos_2sigma_m.powi(2) - 1.0;
        
        // 8. Cálculo de (λ-L) según ecuación 179
        lambda_diff_prev = (e2/2.0) * b0_i.cos() * (
            1.0 + e2/4.0 + e2.powi(2)/8.0 - (e2/8.0 + e2.powi(2)/8.0)*b0_i.sin().powi(2) + (3.0*e2.powi(2)/64.0)*b0_i.sin().powi(4)) * sigma_i
            + ((e2/8.0 + e2.powi(2)/8.0)*b0_i.sin().powi(2) - (e2.powi(2)/16.0)*b0_i.sin().powi(4)) * sigma_i.sin() * cos_2sigma_m
            + (e2.powi(2)/64.0)*b0_i.sin().powi(4) * (2.0*sigma_i).sin() * cos_4sigma_m / 2.0;
        
        
        // 9. Verificación de convergencia
        if iter_count >= max_iter || (lambda_diff_prev - (lambda_prev - delta_lambda)).abs() < threshold {
            if iter_count >= max_iter {
                println!("Advertencia: Máximo de iteraciones alcanzado");
            }
            break;
        }
        
        lambda_prev = delta_lambda + lambda_diff_prev;
    }
    
    // 10. Cálculo final del azimut (ecuación 187)
    let azimuth12 = (delta_lambda.sin() * beta2.cos()).atan2(
        beta2.sin() * beta1.cos() - beta2.cos() * beta1.sin() * delta_lambda.cos()
    ).to_degrees();
    
    // 11. Cálculo final de la distancia geodésica (ecuación 165)
    let two_sigma_m_final = 2.0 * ((beta1.sin() * beta2.sin() / sigma_i.cos()).atan()) - sigma_i;
    let s = b * (
        a_coeff * sigma_i  // Término principal (debería dar ~30,855 m)
        - b_coeff * sigma_i.sin() * two_sigma_m_final.cos()  // Corrección pequeña
        - (c_coeff/2.0) * (2.0*two_sigma_m_final).cos() * (2.0*sigma_i).sin()
    );
    /*let s = b * a_coeff * (sigma_i 
        - b_coeff * sigma_i.sin() * two_sigma_m_final.cos()
        - (c_coeff/2.0) * (2.0*sigma_i).sin() * (2.0*two_sigma_m_final).cos()
    );*/
    let termino_B = b * b_coeff * sigma_i.sin() * two_sigma_m_final.cos();
    let termino_C = b * (c_coeff/2.0) * (2.0*two_sigma_m_final).cos() * (2.0*sigma_i).sin();
    println!("Término s: {} m", s);
    println!("Término B: {} m", termino_B);
    println!("Término C: {} m", termino_C);
    println!("Término principal (b*A*σ): {} m", b * a_coeff * sigma_i);
    println!("b (semieje menor): {} metros", b);
    println!("sigmaii (σ): {} radianes", sigma_i);
    println!("aMayor (A): {}", a_coeff);
    println!("bMayor (B): {}", b_coeff);
    println!("cMayor (C): {}", c_coeff);
        
    (azimuth12, s)
}
fn computeBesselInverso(
    p_geo1: &PointCurvilinear,
    p_geo2: &PointCurvilinear,
    a: f64,
    b: f64
) -> (f64, f64) {
    // 1. Conversión a radianes
    let phi1 = p_geo1.fi.to_radians();
    let phi2 = p_geo2.fi.to_radians();
    let lam1 = p_geo1.la.to_radians();
    let lam2 = p_geo2.la.to_radians();
    
    // 2. Parámetros elipsoidales
    let e2 = (a.powi(2) - b.powi(2)) / a.powi(2); // Excentricidad al cuadrado
    let ep2 = (a.powi(2) - b.powi(2)) / b.powi(2); // Segunda excentricidad al cuadrado
    
    // 3. Latitudes reducidas
    let beta1 = ((1.0 - e2).sqrt() * phi1.tan()).atan();
    let beta2 = ((1.0 - e2).sqrt() * phi2.tan()).atan();
    
    // 4. Diferencia de longitud
    let l_mayus = (lam2 - lam1).abs();
    
    // 5. Cálculo inicial de sigma
    let cos_sigma = beta1.sin() * beta2.sin() + beta1.cos() * beta2.cos() * l_mayus.cos();
    let mut sigma = cos_sigma.acos();
    
    // 6. Azimut inicial
    let alpha = ((beta1.cos() * beta2.cos() * l_mayus.sin()) / sigma.sin()).asin();
    
    // 7. Coeficientes constantes
    let b0 = (alpha.sin() * beta1.cos()).acos();
    let k2 = ep2 * alpha.cos().powi(2);
    
    let a_coeff = 1.0 + (k2/4.0) - (3.0*k2.powi(2)/64.0);
    let b_coeff = (k2/4.0) - (k2.powi(2)/16.0);
    let c_coeff = (k2.powi(2)/64.0);
    
    // 8. Iteración para convergencia
    let mut iter_count = 0;
    let max_iter = 15;
    let threshold = 1e-12;
    let mut sigma_prev;
    
    loop {
        iter_count += 1;
        sigma_prev = sigma;
        
        // 9. Cálculo de términos intermedios
        let two_sigma_m = 2.0 * ((beta1.sin() * beta2.sin() / sigma.cos()).atan()) - sigma;
        let cos_2sigma_m = two_sigma_m.cos();
        
        // 10. Actualización de sigma
        sigma = l_mayus / (a_coeff * b * (1.0 
            - b_coeff * sigma.sin() * cos_2sigma_m 
            - (c_coeff/2.0) * (2.0*sigma).sin() * (2.0*cos_2sigma_m).cos()
        ));
        
        // 11. Verificación de convergencia
        if (sigma - sigma_prev).abs() < threshold || iter_count >= max_iter {
            break;
        }
    }
    
    // 12. Cálculo final de la distancia geodésica
    let two_sigma_m_final = 2.0 * ((beta1.sin() * beta2.sin() / sigma.cos()).atan()) - sigma;
    let s = b * a_coeff * (sigma 
        - b_coeff * sigma.sin() * two_sigma_m_final.cos()
        - (c_coeff/2.0) * (2.0*sigma).sin() * (2.0*two_sigma_m_final).cos()
    );
    
    // 13. Cálculo de azimutes
    let azimuth12 = (alpha.sin() / beta1.cos()).asin();
    let azimuth21 = (alpha.sin() / beta2.cos()).asin();

    
    (azimuth12.to_degrees(), s)
}
/*
fn computeBesselInverso(p_geo1: &PointCurvilinear,p_geo2: &PointCurvilinear,a:f64,b:f64)->(f64,f64){

    let phi1=p_geo1.fi* (PI / 180.0);
    let phi2=p_geo2.fi* (PI / 180.0);
    let lam1=p_geo1.la* (PI / 180.0);
    let lam2=p_geo2.la* (PI / 180.0);
    let delta_lambda = (lam1 - lam2).abs();
    let delta_phi = (phi1 - phi2).abs();
    let e2 = 1.0 - (b.powi(2) / a.powi(2));  
    //Calcular las latitudes reducidas a la esfera
    let b1 =((1.0_f64-e2).sqrt() * phi1.tan()).atan();
    let b2 =((1.0_f64-e2).sqrt() * phi2.tan()).atan();
    
    let lMayus = (lam2 -lam1).abs();
    let sigma = b1.sin() * b2.sin()+ (b1.cos() * b2.cos() * lMayus.cos());
    let azimuth12 = (lMayus.sin() * b2.cos()/sigma.sin()).asin();
    let b0 =(azimuth12.sin()*b1.cos()).acos();
    let aprima = 1.0 
        + e2/4.0 
        + e2.powi(2)/8.0 
        - (e2/8.0 + e2.powi(2)/8.0)*b0.sin().powi(2) 
        + (3.0*e2.powi(2)/64.0)*b0.sin().powi(4);

    let bprima = (e2/8.0 + e2.powi(2)/8.0)*b0.sin().powi(2) 
        - (e2.powi(2)/16.0)*b0.sin().powi(4);

    let cprima = (e2.powi(2)/64.0)*b0.sin().powi(4);
    
    let alpha = ((b1.cos() * b2.cos() * lMayus.sin()) / sigma.sin()).asin(); 
    let cos_2sigma_m = sigma.cos() - (2.0*b1.sin()*b2.sin())/alpha.cos().powi(2);
    let cos_4sigma_m = 2.0*cos_2sigma_m.powi(2) - 1.0;
    let mut l_diff_prev = (e2/2.0) * b0.cos() * (
        aprima * sigma
        + bprima * sigma.sin() * cos_2sigma_m
        + (cprima/2.0) * (2.0*sigma).sin() * cos_4sigma_m);
    let mut lam_prev = lMayus + l_diff_prev;
    let mut iter_count = 0;
    let max_iter = 15;
    let threshold = 0.00001;
    let mut sigmaii=0.0;
    let mut alphaii=0.0;
    let mut b0ii = 0.0;

    loop {
        iter_count += 1;
        if iter_count > max_iter {
            println!("Advertencia: No converge después de {} iteraciones", max_iter);
            break; // Límite de iteraciones razonable
        }
        // Calcular sigma_i usando el lambda actual
        let sigmai= b1.sin() * b2.sin()+ (b1.cos() * b2.cos() * lam_prev.cos());
        sigmaii = sigmai;
        let azimuth12i = (lam_prev.sin() * b2.cos()/sigmai.sin()).asin();
        let b0i =(azimuth12i.sin()*b1.cos()).acos();
        b0ii = b0i;
        let aprimai = 1.0 
        + e2/4.0 
        + e2.powi(2)/8.0 
        - (e2/8.0 + e2.powi(2)/8.0)*b0i.sin().powi(2) 
        + (3.0*e2.powi(2)/64.0)*b0i.sin().powi(4);
        let bprimai = (e2/8.0 + e2.powi(2)/8.0)*b0i.sin().powi(2) 
            - (e2.powi(2)/16.0)*b0i.sin().powi(4);
        let cprimai = (e2.powi(2)/64.0)*b0i.sin().powi(4);
        let alphai = ((b1.cos() * b2.cos() * lam_prev.sin()) / sigmai.sin()).asin(); 
        alphaii = alphai;
        let cos_2sigma_mi = sigmai.cos() - (2.0*b1.sin()*b2.sin())/alphai.cos().powi(2);
        let cos_4sigma_mi = 2.0*cos_2sigma_mi.powi(2) - 1.0;
        let l_diff_new = (e2/2.0) * b0i.cos() * (
            aprimai * sigmai
            + bprimai * sigmai.sin() * cos_2sigma_mi
            + (cprimai/2.0) * (2.0*sigmai).sin() * cos_4sigma_mi);
        if (l_diff_new-l_diff_prev).abs() < threshold{
            break
        }
        l_diff_prev = l_diff_new;
        lam_prev = lMayus + l_diff_prev;
        
    }
    //Calcular azimuths
    let azimuth12d = (alphaii.sin() / b1.cos()).asin();
    let azimuth21d = (alphaii.sin() / b2.cos()).asin();
    let sigma1 = (b1.tan() / azimuth12d.cos()).atan();
    let dosSigmam = 2.0*sigma1+sigmaii;
    //Calcular distancia geodesica s
    let ePrima2 =(a.powi(2)-b.powi(2)) / b.powi(2);
    println!("eprima2  es: {}", ePrima2);
    let k2 = ePrima2 * b0ii.sin().powi(2);
    println!("K  es: {}", k2);
    let aMayor = 1.0+(k2/4.0)-(3.0*(k2.powi(2))/64.0);
    println!("A  es: {}", aMayor);
    let bMayor = (k2/4.0)-(k2.powi(2)/16.0);
    println!("B  es: {}", bMayor);
    let cMayor = (k2.powi(2)/64.0);
    println!("C  es: {}", cMayor);    
    let s = b* ((aMayor*sigmaii)-(bMayor* sigmaii.sin() *dosSigmam.cos())-((cMayor/2.0)*(2.0*dosSigmam).cos()* (2.0*sigmaii).sin()));
    (azimuth12d,s)

}*/


fn main() -> Result<(), Box<dyn Error>>{
    //WGS84 parameters
    let a = 6378137.000;
    let b = 6356752.3142;
    let e2 = 1.0 - ((b as f64).powi(2) / (a as f64).powi(2));
    //Leer correnadas de los puntos de un archivo csv que resulta del software EMLID de los GPS
    /*let puntos=leer_puntos_archivo("C:\\Users\\rob_r\\OneDrive\\Documentos\\geodesiaNotas\\geodesia\\geodesiaII\\datosLevantamiento27Marzo\\gps.csv")?;
    let mut pointsCartesianos = Vec::new();
    let mut distanciasPuntos = Vec::new();
    for (i,point) in puntos.iter().enumerate(){
        let npunto = computeN(a,b,&point);
        let vpunto = computeCartesian(&point,npunto,a,b);
        pointsCartesianos.push(vpunto);
       println!("Punto {}: latitud={}, latitud={}, altura elipsoidal={}", i+1,point.fi,point.la,point.h);
    }
    for i in 0..pointsCartesianos.len(){
       let puntoactual= pointsCartesianos[i];
        let siguientePunto = if i == pointsCartesianos.len()-1{
            pointsCartesianos[0]
        } else{
            pointsCartesianos[i+1]
        };
        let resta = puntoactual.subtract(&siguientePunto);
       let dist = resta.module();
        distanciasPuntos.push(dist);
        println!("Punto {}: x={}, y={},  z={}", i+1,puntoactual.x,puntoactual.y,puntoactual.z);
        println!("Distancia entre punto {} y {}: {:.4}", i+1, 
                 if i == pointsCartesianos.len() - 1 { 1 } else { i+2 }, 
                 dist);
    }
    Fin de el segmento que lee las coordenadas de los puntos de un archivo csv y calcula la distancia entre los puntos*/

    
    //Aqui empieza el segmento para calcular distancias geodesicas con Bessel
    //Problema directo con Bessel primer punto, distancia y azimuth al segundo punto conocidos
    /*
    let p1  =read_coordinate("Ingrese las coordenadas del primer punto");
    let azdist = read_azdist("Ingrese la distancia y azimuth al siguiente punto");
    let azimuth=azdist.azimuth12;
    let azimuthRadians = azimuth.to_radians();
    let s= azdist.dist12;
    let n1 = computeN(a,b,&p1);
    //Paso 1 Calcular la latitud reducida Beta1 
    let latB1=latitudB1(&p1,a, b);
    println!("La latitud reducida B1  es: {}", latB1* (180.0/PI ));
    //Paso 2 calcular el azimuth de la geodesica en el ecuador
    let azec = azimuthRadians.sin() * latB1.cos();
    println!("Azimuth de la geodesica en el ecuador  es: {}", azec);
    //Paso 3 calcular sigma a traves de iteraciones
    let sigmaT = iterarSigma(a,b, s, latB1, azimuthRadians);
    //Paso 4 calcular B2 usando 145 y B0 con 143a
    let sig1 = (latB1.tan() / azimuthRadians.cos()).atan();
    let latB0 =(azimuthRadians.sin()*latB1.cos()).acos();
    let latB2 =((sig1+sigmaT).sin()*(latB0).sin()).asin();
    //Paso 5 Calcular la latitud del punto 2 con la ecuacion 128
    let latitud1punto =p1.la* (PI /180.0);
    let fi2= (((latB2).tan())/((1.0_f64-e2).sqrt())).atan();
    println!("La latitud del punto 2 Bessel  es: {}", fi2* (180.0/PI ));
    //Paso 6 calcular la longitud del segundo punto con 179 iterando
    let (lon2pre,alphaecuador) = lonBessel(a,b,sigmaT,azimuthRadians,latB0,latB1,latB2);
    let lon2 = latitud1punto+ lon2pre;
    println!("La longitud del punto 2 Bessel es: {}", lon2* (180.0/PI ));
    let azimuth_inverso = reverse_azimuthBessel(alphaecuador,latB2);
    println!("el azimuth inverso es: {}", azimuth_inverso);
    let p2 =compute_p2_Puissant(&p1, a, b, s, azimuth);
    println!("La latitud del punto 2  es: {}", p2.fi);
    println!("La longitud del punto 2  es: {}", p2.la);*/
    //Aqui termina el problema directo con Bessel


     //Aqui empieza el problema inverso con Bessel
     //Pedir al usuario dos coordenadas 
    let p1  =read_coordinate("Ingrese las coordenadas del primer punto");
    let p2  =read_coordinate("Ingrese las coordenadas del segundo punto");
    let (resultadoAzimuth,resultadodistancia)=computePuissantInverso(&p1,&p2, a, b);
    println!("Azimuth con Puissant  es: {}", resultadoAzimuth);
    println!("Distancia Geodesica con Puissant  es: {}", resultadodistancia);
    let (resultadoAzimuthBessel,resultadodistanciaBessel)= compute_bessel_inverse(&p1,&p2, a, b);
    println!("Azimuth con Bessel  es: {}", resultadoAzimuthBessel);
    println!("Distancia Geodesica con Bessel  es: {}", resultadodistanciaBessel);



    // Pedir al usuario los dos vectores
    //let v1 = read_vector("Ingrese el primer vector:");
    //let v2 = read_vector("Ingrese el segundo vector:");

    //Pedir al usuario dos coordenadas 
    //let p1  =read_coordinate("Ingrese las coordenadas del primer punto");
    // let p2  =read_coordinate("Ingrese las coordenadas del segundo punto");


    //let (resultadoAzimuth,resultadodistancia)=computePuissantInverso(&p1,&p2, a, b);

    //problema directo punto , distancia y azimuth
    //let azdist = read_azdist("Ingrese la distancia y azimuth al siguiente punto");
    //let azimuth=azdist.azimuth12;
    //let s= azdist.dist12;
    //let p2  =read_coordinate("Ingrese las coordenadas del segundo punto");
    //Calcular los vectores Normal de cada punto 
    //let n1 = computeN(a,b,&p1);
    //let p2 =compute_p2_Puissant(&p1, a, b, s, azimuth);
    //println!("La latitud del punto 2  es: {}", p2.fi);
    //println!("La longitud del punto 2  es: {}", p2.la);
    //Calcular la distancia entre ambos puntos para comprobar aproximacion, usando los arcos de meridiano y de paralelo
    //let arcoMer =arcoMeridiano(&p1, &p2, a, b);
    //let arcoMer =arcoMeridianoPromedio(&p1, &p2, a, b);
    //let arcoPar =arcoParalelo(&p1, &p2, a, b);
    //let distanciaPuntos = calcularDistanciaArcos(arcoMer,arcoPar);
    //println!("El arco de meridiano entre ambos puntos es de : {}", arcoMer);
    //println!("El arco de paralelo entre ambos puntos es de {}", arcoPar);
    //println!("La distancia entre ambos puntos usando los arcos de meridiano y paralelo es de  {}", distanciaPuntos);
    //println!("La diferencia entre la distancia original y la calculada es de {}", s-distanciaPuntos);
    //let n2 = computeN(a,b,&p2);
    //let v1 = computeCartesian(&p1,n1,a,b);
    //let v2 = computeCartesian(&p2,n2,a,b);
    //let resultado_resta = v1.subtract(&v2);
    //let modulo = resultado_resta.module();
    //println!("La distancia entre ambos puntos usando la distancia euclideana es de  {}", modulo);
    //println!("La diferencia entre la distancia original y la geometrica es de {}", s-modulo);



    //let n2 = computeN(a,b,&p2);
    //Obtener las coordenadas cartesianas de cada punto
    //let v1 = computeCartesian(&p1,n1,a,b);
    //let v2 = computeCartesian(&p2,n2,a,b);
    //let resultado_resta = v1.subtract(&v2);
    //let modulo = resultado_resta.module();
    //println!("La primer vertical n1 es: {}", n1);
    //println!("La segunda vertical n1 es: {}", n2);
    //println!("Las coordenadas cartesianas del primer punto son: {:?}", v1);
    //println!("Las coordenadas cartesianas del segundo punto son: {:?}", v2);
    //println!("El vector formado de ambos puntos es: {:?}", resultado_resta);
    //println!("El módulo del vector resultante es: {}", modulo);
    //println!("Calculando las coordenadas geodesicas de vuelta");
    //let longitud1 = computeLongitude(&v1);
    //println!("longitud1: {}", longitud1);
    
    //let longitud2 =computeLongitude(&v2);
    //println!("longitud2: {}", longitud2);
    //let (h11, phi11, n11) = compute_iterative(&v1, a, b);
    //println!("h: {}", h11);
    //println!("phi: {} radianes", phi11);
    //println!("phi: {} grados", phi11 * 180.0 / PI);
    //println!("n: {}", n11);
    //let (h22, phi22, n22) = compute_iterative(&v2, a, b);
    //println!("h: {}", h22);
    //println!("phi: {} radianes", phi22);
    //println!("phi: {} grados", phi22 * 180.0 / PI);
    //println!("n: {}", n22);
    Ok(())



}