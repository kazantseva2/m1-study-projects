package org.vmk.dep508.io.rest;

import org.apache.log4j.Logger;

import java.io.*;
import java.net.Inet4Address;
import java.net.InetAddress;
import java.net.Socket;

/*
POST /books/1491927925 HTTP/1.1
cache-control: no-cache
content-type: application/x-www-form-urlencoded

{"isbn":"1491927925","title":"MapReduce Design Patterns 2nd edition","authors":["Donald Miner", "Adam Shook"],"publisher":"O'Reilly Media","annotation":"<p>Until now, design patterns for the MapReduce framework have been scattered among various research papers, blogs, and books. This handy guide brings together a unique collection of valuable MapReduce patterns that will save you time and effort regardless of the domain, language, or development framework you’re using.</p><p>Each pattern is explained in context, with pitfalls and caveats clearly identified to help you avoid common design mistakes when modeling your big data architecture. This book also provides a complete overview of MapReduce that explains its origins and implementations, and why design patterns are so important. All code examples are written for Hadoop.</p>"}



{"isbn":"1491901632","title":"Hadoop: The Definitive Guide","authors":["Tom White"],"publisher":"O'Reilly Media","annotation":"<p>Get ready to unlock the power of your data. With the fourth edition of this comprehensive guide, you’ll learn how to build and maintain reliable, scalable, distributed systems with Apache Hadoop. This book is ideal for programmers looking to analyze datasets of any size, and for administrators who want to set up and run Hadoop clusters. Using Hadoop 2 exclusively, author Tom White presents new chapters on YARN and several Hadoop-related projects such as Parquet, Flume, Crunch, and Spark. You’ll learn about recent changes to Hadoop, and explore new case studies on Hadoop’s role in healthcare systems and genomics data processing.</p>"}

{"isbn":"1449358624","title":"Learning Spark: Lightning-Fast Big Data Analysis","authors":["Holden Karau", "Andy Konwinski", "Patrick Wendell", "Matei Zaharia"],"publisher":"O'Reilly Media","annotation":"<p>Data in all domains is getting bigger. How can you work with it efficiently? Recently updated for Spark 1.3, this book introduces Apache Spark, the open source cluster computing system that makes data analytics fast to write and fast to run. With Spark, you can tackle big datasets quickly through simple APIs in Python, Java, and Scala. This edition includes new information on Spark SQL, Spark Streaming, setup, and Maven coordinates.</p>"}


 */

public class Client {

    private static Logger log = Logger.getLogger(Client.class);

    public static void main(String[] args) throws IOException {
        saveBooks();
        deleteBook();
    }

    private static void deleteBook() throws IOException {
        InetAddress address = Inet4Address.getByName("localhost");
        int port = 8080;

        Socket socket = new Socket(address, port);
        BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(socket.getOutputStream()));

        writer.write("DELETE /books/1449358624 HTTP/1.1");
        writer.newLine();
        writer.write("cache-control: no-cache");
        writer.newLine();
        writer.write("content-type: application/x-www-form-urlencoded");
        writer.newLine();
        writer.newLine();

        writer.flush();

        BufferedReader reader = new BufferedReader(new InputStreamReader(socket.getInputStream()));
        reader.lines().forEach(log::info);
    }

    private static void saveBooks() throws IOException {
        InetAddress address = Inet4Address.getByName("localhost");
        int port = 8080;

        Socket socket = new Socket(address, port);
        BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(socket.getOutputStream()));

        writer.write("POST /books/1491927925 HTTP/1.1");
        writer.newLine();
        writer.write("cache-control: no-cache");
        writer.newLine();
        writer.write("content-type: application/x-www-form-urlencoded");
        writer.newLine();
        writer.newLine();
        writer.write("{\"isbn\":\"1491927925\",\"title\":\"MapReduce Design Patterns 2nd edition\",\"authors\":[\"Donald Miner\", \"Adam Shook\"],\"publisher\":\"O'Reilly Media\",\"annotation\":\"<p>Until now, design patterns for the MapReduce framework have been scattered among various research papers, blogs, and books. This handy guide brings together a unique collection of valuable MapReduce patterns that will save you time and effort regardless of the domain, language, or development framework you’re using.</p><p>Each pattern is explained in context, with pitfalls and caveats clearly identified to help you avoid common design mistakes when modeling your big data architecture. This book also provides a complete overview of MapReduce that explains its origins and implementations, and why design patterns are so important. All code examples are written for Hadoop.</p>\"}");
        writer.newLine();
        writer.flush();

        writer.write("POST /books/1491901632 HTTP/1.1");
        writer.newLine();
        writer.write("cache-control: no-cache");
        writer.newLine();
        writer.write("content-type: application/x-www-form-urlencoded");
        writer.newLine();
        writer.newLine();
        writer.write("{\"isbn\":\"1491901632\",\"title\":\"Hadoop: The Definitive Guide\",\"authors\":[\"Tom White\"],\"publisher\":\"O'Reilly Media\",\"annotation\":\"<p>Get ready to unlock the power of your data. With the fourth edition of this comprehensive guide, you’ll learn how to build and maintain reliable, scalable, distributed systems with Apache Hadoop. This book is ideal for programmers looking to analyze datasets of any size, and for administrators who want to set up and run Hadoop clusters. Using Hadoop 2 exclusively, author Tom White presents new chapters on YARN and several Hadoop-related projects such as Parquet, Flume, Crunch, and Spark. You’ll learn about recent changes to Hadoop, and explore new case studies on Hadoop’s role in healthcare systems and genomics data processing.</p>\"}");
        writer.newLine();
        writer.flush();

        writer.write("POST /books/1449358624 HTTP/1.1");
        writer.newLine();
        writer.write("cache-control: no-cache");
        writer.newLine();
        writer.write("content-type: application/x-www-form-urlencoded");
        writer.newLine();
        writer.newLine();
        writer.write("{\"isbn\":\"1449358624\",\"title\":\"Learning Spark: Lightning-Fast Big Data Analysis\",\"authors\":[\"Holden Karau\", \"Andy Konwinski\", \"Patrick Wendell\", \"Matei Zaharia\"],\"publisher\":\"O'Reilly Media\",\"annotation\":\"<p>Data in all domains is getting bigger. How can you work with it efficiently? Recently updated for Spark 1.3, this book introduces Apache Spark, the open source cluster computing system that makes data analytics fast to write and fast to run. With Spark, you can tackle big datasets quickly through simple APIs in Python, Java, and Scala. This edition includes new information on Spark SQL, Spark Streaming, setup, and Maven coordinates.</p>\"}");
        writer.newLine();
        writer.flush();

        BufferedReader reader = new BufferedReader(new InputStreamReader(socket.getInputStream()));
        reader.lines().forEach(log::info);
    }
}
