package org.vmk.dep508.library;

import org.h2.tools.RunScript;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.vmk.dep508.library.Library;
import org.vmk.dep508.library.LibraryImpl;

import java.io.FileReader;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.Statement;
import java.util.List;

import static org.junit.Assert.*;

public class LibraryTest {
    private Library library;
    private Connection connection;

    @Before
    public void setUp() throws Exception {
        library = new LibraryImpl("jdbc:h2:mem:library", "", "");
        connection = DriverManager.getConnection("jdbc:h2:mem:library");
        RunScript.execute(connection, new FileReader("src/main/java/org/vmk/dep508/library/tables.sql"));
    }

    @After
    public void tearDown() throws Exception {
        try(Statement stmt = connection.createStatement();) {
            String tableSql = "drop table books; drop table abonents;";
            stmt.execute(tableSql);
        }
        connection.close();
    }

    @Test
    public void addNewBook() throws Exception {
        Book newBook1 = new Book(1, "Джейн Эйр");
        Book newBook2 = new Book(2, "Архипелаг ГУЛАГ");
        library.addNewBook(newBook1);
        library.addNewBook(newBook2);

        List<Book> availableBooks = library.findAvailableBooks();
        assertEquals(2, availableBooks.size());
        assertTrue(availableBooks.contains(newBook1));
        assertTrue(availableBooks.contains(newBook2));
    }

    @Test
    public void addAbonent() throws Exception {
        Student student = new Student(1, "Varvara");
        library.addAbonent(student);

        List<Student> students = library.getAllStudents();
        assertEquals(1, students.size());
        assertEquals(student, students.get(0));
    }

    @Test
    public void borrowBook() throws Exception {
        Book book1 = new Book(1, "Джейн Эйр");
        Book book2 = new Book(2, "Архипелаг ГУЛАГ");
        Student student1 = new Student(1, "Varvara");
        Student student2 = new Student(2, "Tatyana");

        library.addNewBook(book1);
        library.borrowBook(book1, student1);

        List<Book> availableBooks1 = library.findAvailableBooks();
        assertFalse(availableBooks1.contains(book1));

        List<Student> students = library.getAllStudents();
        assertTrue(students.contains(student1));

        library.addNewBook(book2);
        library.addAbonent(student2);
        library.borrowBook(book2, student2);

        List<Book> availableBooks2 = library.findAvailableBooks();
        assertFalse(availableBooks2.contains(book2));
    }

    @Test
    public void returnBook() throws Exception {
        Book book = new Book(1, "Джейн Эйр");
        Student student = new Student(1, "Varvara");

        library.addNewBook(book);
        library.borrowBook(book, student);
        library.returnBook(book, student);

        List<Book> availableBooks = library.findAvailableBooks();
        assertTrue(availableBooks.contains(book));
    }

    @Test
    public void findAvailableBooks() throws Exception {
        Book book1 = new Book(1, "Джейн Эйр");
        Book book2 = new Book(2, "Архипелаг ГУЛАГ");
        library.addNewBook(book1);
        library.addNewBook(book2);

        List<Book> availableBooks = library.findAvailableBooks();
        assertEquals(2, availableBooks.size());
        assertTrue(availableBooks.contains(book1));
        assertTrue(availableBooks.contains(book2));
    }

    @Test
    public void getAllStudents() throws Exception {
        Student student1 = new Student(1, "Varvara");
        Student student2 = new Student(2, "Tatyana");

        library.addAbonent(student1);
        library.addAbonent(student2);

        List<Student> students = library.getAllStudents();
        assertEquals(2, students.size());
        assertTrue(students.contains(student1));
        assertTrue(students.contains(student2));
    }

}