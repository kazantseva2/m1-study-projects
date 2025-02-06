package org.vmk.dep508.stream.lenses;


import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

public class App {
    public static void main(String[] args) throws IOException {
        Files
                .lines(Paths.get("lenses.data"))
                .map(LenseData::parse)
                .forEach(System.out::println);
    }
}
