fn main() {
    cc::Build::new()
        .file("extern/hke-sdg/keplerh.c")
        .file("extern/hke-sdg/polinomios25.c")
        .file("extern/hke-sdg/polinomios25Q.c")
        .opt_level(3)
        .compile("hke-sdg");
}
