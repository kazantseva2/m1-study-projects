package org.vmk.dep508.io;

import com.fasterxml.jackson.databind.ObjectMapper;
import org.apache.log4j.Logger;
import org.vmk.dep508.io.rest.Book;

import java.io.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class SerializationExample {

    private static Logger log = Logger.getLogger(SerializationExample.class);

    public static void main(String[] args) throws IOException, ClassNotFoundException {

        Book b = new Book(
                "1491927925",
                "MapReduce Design Patterns 2nd edition",
                new String[] {"Donald Miner"},
                "O'Reilly Media",
                "<p>Until now, design patterns for the MapReduce framework have been scattered among various research papers, blogs, and books. This handy guide brings together a unique collection of valuable MapReduce patterns that will save you time and effort regardless of the domain, language, or development framework you’re using.</p>" +
                        "<p>Each pattern is explained in context, with pitfalls and caveats clearly identified to help you avoid common design mistakes when modeling your big data architecture. This book also provides a complete overview of MapReduce that explains its origins and implementations, and why design patterns are so important. All code examples are written for Hadoop.</p>");

        //serialize object and save to file
        try (ObjectOutputStream output =
                     new ObjectOutputStream(
                             new FileOutputStream("1491927925.data"))){
            output.writeObject(b);
            output.flush();
        }

        //serialize object, gzip and save to file
        try (ObjectOutputStream output =
                     new ObjectOutputStream(
                             new GZIPOutputStream(
                                     new FileOutputStream("1491927925.gzip")))){
            output.writeObject(b);
            output.flush();
        }

        //serialize object into JSON and save to file
        try (OutputStream output = new FileOutputStream("1491927925.json")) {
            ObjectMapper mapper = new ObjectMapper();
            mapper.writeValue(output, b);
            output.flush();
        }

        try (ObjectInputStream input =
                     new ObjectInputStream(
                             new GZIPInputStream(
                                     new FileInputStream("1491927925.gzip")))) {
            Book unserialized = (Book)input.readObject();
            log.info(unserialized);
        }
    }
}
