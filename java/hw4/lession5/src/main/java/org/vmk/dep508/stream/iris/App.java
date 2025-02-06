package org.vmk.dep508.stream.iris;

import org.vmk.dep508.stream.employees.Employee;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;


public class App {
    public static void main(String[] args) throws IOException {
        App a = new App();
        a.test();
    }

    public void test() throws IOException {

        //load data from file iris.data
        List<Iris> irisList = Files
                .lines(Paths.get("iris.data"))
                .map(Iris::parse)
                .collect(Collectors.toList());
        IrisDataSetHelper helper = new IrisDataSetHelper(irisList);

        //get average sepal width
        Double avgSepalLength = helper.getAverage(Iris::getSepalWidth);
        System.out.println(avgSepalLength);

        //get average petal square - petal width multiplied on petal length
        Double avgPetalLength = helper.getAverage(iris -> iris.getPetalWidth() * iris.getPetalLength());
        System.out.println(avgPetalLength);

        //get average petal square for flowers with sepal width > 4
        Double avgPetalSquare = helper.getAverageWithFilter(
                iris -> iris.getSepalWidth() > 4,
                iris -> iris.getPetalWidth() * iris.getPetalLength());
        System.out.println(avgPetalSquare);

        //get flowers grouped by Petal size (Petal.SMALL, etc.)
        Map groupsByPetalSize = helper.groupBy(Iris::classifyByPatel);
        System.out.println(groupsByPetalSize);

        //get max sepal width for flowers grouped by species
        Map maxSepalWidthForGroupsBySpecies = helper.maxFromGroupedBy(Iris::getSpecies, Iris::getSepalWidth);
        System.out.println(maxSepalWidthForGroupsBySpecies);
    }
}

